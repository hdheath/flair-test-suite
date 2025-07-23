# src/flair_test_suite/qc/region_qc.py
# -------------------------------------------------
# QC collector for the **slice** stage.
#
# This module summarises every per‑region folder created by `SliceStage`.
# For each region it reports
#   • number of alignments in the sliced BAM
#   • number of lines in the sliced BED
#   • (optional) line‑counts for any extra files that were requested
#   • basic gene / transcript statistics parsed directly from the small
#     region‑specific GTF (no gffutils or SQLite overhead!)
#
# Two artefacts are written into the slice stage directory:
#   1) <stage_dir>/slice_region_details.tsv   – one row *per region*
#   2) <stage_dir>/slice_region_qc.tsv        – run‑level summary (row‑count,
#                                               QC runtime, slice runtime)
#
# If an individual slice ends up completely empty (no BAM/BED/GTF data)
# a `UserWarning` is emitted so the caller can decide whether to flag that.

from __future__ import annotations
from pathlib import Path
import subprocess
import csv
import time
from typing import List, Dict, Any
import warnings

from . import register, write_metrics  # registry + helper from qc/__init__.py
from .qc_utils import count_lines      # shared cheap wc -l helper

__all__ = ["collect"]

# ---------------------------------------------------------------------------
# Constants & helper functions
# ---------------------------------------------------------------------------

# Map file‑name suffix ➜ column prefix in details table
OPTIONAL_EXTS: dict[str, str] = {
    "tab":       "junctions",   # STAR SJ.out.tab
    "exp5.bed":  "exp5_bed",
    "exp3.bed":  "exp3_bed",
    "ref5.bed":  "ref5_bed",
    "ref3.bed":  "ref3_bed",
}


def _bam_count(bam: Path) -> int:  # pragma: no cover (tiny helper)
    """Return `samtools view -c` count for **primary** alignments.

    We do not trap non‑zero exit codes – if samtools is missing or the BAM is
    corrupt we *want* the pipeline to fail loudly instead of silently masking
    the problem.
    """
    if not bam.exists():
        return 0
    return int(
        subprocess.check_output(["samtools", "view", "-c", str(bam)], text=True).strip()
    )


def _parse_attrs(attr_col: str) -> Dict[str, str]:
    """Parse the 9th GTF column (attributes) into a dict."""
    attrs: Dict[str, str] = {}
    for chunk in attr_col.rstrip(";").split(";"):
        chunk = chunk.strip()
        if chunk:
            if " " not in chunk:
                continue  # malformed piece – skip with no warning
            k, v = chunk.split(" ", 1)
            attrs[k] = v.strip('"')
    return attrs


def _gene_tx_summary(gtf: Path) -> Dict[str, Any]:
    """Compute lightweight per‑slice transcript / gene metrics.

    Because each region GTF is tiny (a few hundred lines at most), we can scan
    it once without building any external database.
    """
    if not gtf.exists() or gtf.stat().st_size == 0:
        # Empty slice – return zeros so caller can still write a row
        return {
            "gene_count": 0,
            "transcript_count": 0,
            "mean_tx_len": 0,
            "median_tx_len": 0,
            "mean_exons": 0,
            "median_exons": 0,
            "mean_iso_per_gene": 0,
            "median_iso_per_gene": 0,
        }

    # --- first pass – gather per‑transcript info ---------------------------
    tx_lengths: Dict[str, int] = {}
    exon_counts: Dict[str, int] = {}
    gene_to_tx: Dict[str, List[str]] = {}

    with open(gtf) as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                # malformed line – warn once, then continue
                warnings.warn(f"Malformed GTF line skipped in {gtf.name}: {raw[:50]}…", UserWarning)
                continue
            feature, start_s, end_s, attrs_s = cols[2], cols[3], cols[4], cols[8]
            try:
                start_i, end_i = int(start_s), int(end_s)
            except ValueError:
                warnings.warn(f"Non‑integer coords in {gtf.name}: {start_s}/{end_s}", UserWarning)
                continue
            at = _parse_attrs(attrs_s)
            tid = at.get("transcript_id")
            gid = at.get("gene_id")

            if feature == "transcript" and tid:
                tx_lengths[tid] = end_i - start_i + 1
                if gid:
                    gene_to_tx.setdefault(gid, []).append(tid)
            elif feature == "exon" and tid:
                exon_counts[tid] = exon_counts.get(tid, 0) + 1

    # no transcripts → no meaningful stats
    if not tx_lengths:
        return {
            "gene_count": 0,
            "transcript_count": 0,
            "mean_tx_len": 0,
            "median_tx_len": 0,
            "mean_exons": 0,
            "median_exons": 0,
            "mean_iso_per_gene": 0,
            "median_iso_per_gene": 0,
        }

    # --- helper lambdas ----------------------------------------------------
    def _mean(vals: List[int]) -> float:
        return round(sum(vals) / len(vals), 2) if vals else 0.0

    def _median(vals: List[int]) -> float:
        if not vals:
            return 0.0
        vals = sorted(vals)
        mid = len(vals) // 2
        return float(vals[mid]) if len(vals) % 2 else round((vals[mid - 1] + vals[mid]) / 2, 2)

    # --- derive vectors ----------------------------------------------------
    tx_len_vec = list(tx_lengths.values())
    exon_vec   = [exon_counts.get(tid, 0) for tid in tx_lengths]
    iso_vec    = [len(txs)           for txs in gene_to_tx.values()]

    return {
        "gene_count":           len(gene_to_tx),
        "transcript_count":    len(tx_lengths),
        "mean_tx_len":         _mean(tx_len_vec),
        "median_tx_len":       _median(tx_len_vec),
        "mean_exons":          _mean(exon_vec),
        "median_exons":        _median(exon_vec),
        "mean_iso_per_gene":   _mean(iso_vec),
        "median_iso_per_gene": _median(iso_vec),
    }

# ---------------------------------------------------------------------------
# Collector entry‑point – registered for the "slice" stage
# ---------------------------------------------------------------------------

@register("slice")
def collect(manifest: Path, out_dir: Path, runtime_sec: float | None = None):
    """Aggregate QC for *all* regions produced by one `SliceStage`.

    Parameters
    ----------
    manifest     path to `slices_manifest.tsv` (written by the stage)
    out_dir      the slice stage directory – side‑cars are written here
    runtime_sec  wall‑time (sec) spent in `SliceStage.run()` – forwarded so
                 we can surface it alongside our own QC runtime
    """
    qc_start = time.time()

    # ------------------------ gather per‑region rows ----------------------
    rows: List[Dict[str, Any]] = []
    if not manifest.exists():
        raise FileNotFoundError(f"slice manifest not found: {manifest}")

    with open(manifest) as fh:
        header = next(fh).rstrip("\n").split("\t")
        if header[:4] != ["chrom", "start", "end", "region_dir"]:
            warnings.warn("Unrecognised slice manifest header – columns may be mis‑aligned", UserWarning)
        for lineno, raw in enumerate(fh, start=2):
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 4:
                warnings.warn(f"Manifest line {lineno} malformed – skipping", UserWarning)
                continue
            chrom, start_s, end_s, region_dir_s = parts[:4]
            region_dir = Path(region_dir_s)
            tag = f"{chrom}_{start_s}_{end_s}"

            # --- mandatory files -----------------------------------------
            bam = region_dir / f"{tag}.bam"
            bed = region_dir / f"{tag}.bed"
            gtf = region_dir / f"{tag}.gtf"

            rec: Dict[str, Any] = {
                "region":     tag,
                "bam_reads":  _bam_count(bam),
                "bed_lines":  count_lines(bed) if bed.exists() else 0,
            }
            rec.update(_gene_tx_summary(gtf))

            # --- optional artefacts --------------------------------------
            for suffix, col_prefix in OPTIONAL_EXTS.items():
                f = next(region_dir.glob(f"*{suffix}"), None)
                rec[f"{col_prefix}_lines"] = count_lines(f) if f and f.exists() else 0

            rows.append(rec)

            # warn if slice appears empty
            if rec["bam_reads"] == 0 and rec["bed_lines"] == 0 and rec["gene_count"] == 0:
                warnings.warn(
                    f"Slice region {tag} appears empty (no BAM/BED/GTF records)",
                    UserWarning,
                )

    # ------------------------ write details table ------------------------
    details_tsv = out_dir / "slice_region_details.tsv"
    if rows:
        with open(details_tsv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    else:
        warnings.warn("No region rows collected – details table not written", UserWarning)

    # ------------------------ write run‑level side‑car -------------------
    metrics = {
        "region_count":   len(rows),
        "qc_runtime_sec": round(time.time() - qc_start, 2),
    }
    if runtime_sec is not None:
        metrics["slice_runtime_sec"] = round(runtime_sec, 2)

    write_metrics(out_dir, "slice_region", metrics)
    return metrics
