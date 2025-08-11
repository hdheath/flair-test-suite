from __future__ import annotations
from pathlib import Path
from typing import Optional, Dict
import time, json, logging

from . import register, write_metrics
from .qc_utils import (
    count_lines,
    percent,
    count_unique_junctions,
    count_splice_junction_motifs,
)
from ..lib.signature import qc_sidecar_path

__all__ = ["collect", "run_qc"]

def _find_align_stage_dir(out_dir: Path, align_sig: str) -> Path:
    cur = out_dir
    for _ in range(6):
        cand = cur / "align" / align_sig
        if cand.exists():
            return cand
        cur = cur.parent
    raise FileNotFoundError(f"Could not locate align stage dir from {out_dir} (sig={align_sig})")

def _bed_is_empty(p: Path) -> bool:
    try:
        if (not p.exists()) or (p.stat().st_size == 0):
            return True
        return count_lines(p) == 0
    except Exception as e:
        logging.warning(f"[correct-qc] Failed to probe BED {p} ({e}); treating as empty.")
        return True

@register("correct")
def collect(
    bed: Path,
    out_dir: Path,
    n_input_reads: Optional[int],
    align_sig: Optional[str],
    genome_fa: str,
    runtime_sec: float | None = None,
    regionalized: bool = False,
) -> Dict:
    """
    Regionalized=True  -> NO align comparison. Only 'after' motifs from the corrected BED.
    Regionalized=False -> Full QC with before/after/diff and removal metrics.
    """
    logging.info(f"[correct-qc] bed={bed} out_dir={out_dir} regionalized={regionalized}")
    t0 = time.time()

    # Always skip empties early
    if _bed_is_empty(bed):
        metrics = {"skipped_empty_bed": True, "n_corrected_reads": 0, "qc_runtime_sec": 0.0}
        write_metrics(out_dir, "correct", metrics)
        return metrics

    # ────────────── Regionalized: after-only, no align usage ──────────────
    if regionalized:
        try:
            raw_after = count_splice_junction_motifs(bed_path=bed, fasta_path=Path(genome_fa), max_workers=4)
            motif_after = {f"{k[0]}:{k[1]}": v for k, v in raw_after.items()}
        except Exception as e:
            logging.error(f"[correct-qc] Motif counting failed (regionalized): {e}")
            motif_after = {}

        # One JSON per region dir, with AFTER only
        out_json = Path(out_dir) / "correct_splice_site_motifs.json"
        try:
            with open(out_json, "w") as fh:
                json.dump({"after": motif_after}, fh, indent=2)
        except Exception as e:
            logging.error(f"[correct-qc] Failed writing regional motif JSON: {e}")

        metrics = {
            "n_corrected_reads": count_lines(bed),
            "qc_runtime_sec": round(time.time() - t0, 2),
        }
        write_metrics(out_dir, "correct", metrics)
        return metrics

    # ────────────── Non-regionalized: full QC vs align ──────────────
    # 1) counts
    try:
        n_corr = count_lines(bed)
    except Exception as e:
        logging.error(f"[correct-qc] Counting lines failed: {e}")
        n_corr = 0

    if n_input_reads is None:
        n_removed = None
        removed_pct = None
        logging.warning("[correct-qc] n_input_reads is None for non-regionalized run; removal % omitted.")
    else:
        n_removed = n_input_reads - n_corr
        removed_pct = percent(n_removed, n_input_reads)

    # 2) unique junctions after
    try:
        uniq_after = count_unique_junctions(bed)
    except Exception as e:
        logging.error(f"[correct-qc] unique_junctions(after) failed: {e}")
        uniq_after = 0

    # 3) load align QC + motifs (BEFORE)
    uniq_before = None
    motifs_before: dict = {}
    try:
        if not align_sig:
            raise RuntimeError("align_sig is required for non-regionalized QC")
        align_stage_dir = _find_align_stage_dir(out_dir, align_sig)
        qc_tsv = qc_sidecar_path(align_stage_dir, "align")
        if align_stage_dir.is_dir() and qc_tsv.exists():
            with open(qc_tsv) as fh:
                next(fh, None)
                for line in fh:
                    k, v = line.rstrip("\n").split("\t")
                    if k == "unique_junctions":
                        uniq_before = int(v); break
        motif_json_path = align_stage_dir / "splice_site_motifs.json"
        if motif_json_path.exists():
            with open(motif_json_path) as fh:
                motifs_before = json.load(fh)
    except Exception as e:
        logging.error(f"[correct-qc] Failed loading BEFORE from align: {e}")

    # 4) motifs AFTER
    try:
        raw_after = count_splice_junction_motifs(bed_path=bed, fasta_path=Path(genome_fa), max_workers=4)
        motifs_after = {f"{k[0]}:{k[1]}": v for k, v in raw_after.items()}
    except Exception as e:
        logging.error(f"[correct-qc] Motif counting failed (after): {e}")
        motifs_after = {}

    # 5) diff + JSON
    motif_diff = {k: motifs_after.get(k, 0) - motifs_before.get(k, 0)
                  for k in set(motifs_after) | set(motifs_before)}
    out_json = Path(out_dir) / "correct_splice_site_motifs.json"
    try:
        with open(out_json, "w") as fh:
            json.dump({"before": motifs_before, "after": motifs_after, "diff": motif_diff}, fh, indent=2)
    except Exception as e:
        logging.error(f"[correct-qc] Failed writing correct_splice_site_motifs.json: {e}")

    metrics = {
        "n_input_reads": n_input_reads,
        "n_corrected_reads": n_corr,
        "n_removed_reads": n_removed,
        "removed_pct": removed_pct,
        "unique_junc_before": uniq_before,
        "unique_junc_after": uniq_after,
        "correct_runtime_sec": round(runtime_sec, 2) if runtime_sec is not None else None,
        "qc_runtime_sec": round(time.time() - t0, 2),
    }
    write_metrics(out_dir, "correct", metrics)
    return metrics


def run_qc(
    bed_files: list[tuple[Path, str]],
    stage_dir: Path,
    n_input_reads: Optional[int],
    align_sig: Optional[str],
    genome_fa: str,
    runtime_sec: float | None,
    is_regionalized: bool,
) -> dict:
    """
    Iterate regions, call collect() with explicit 'regionalized' flag.
    - Regionalized: QC outputs go under <stage_dir>/<region_tag>/
    - Non-regionalized: QC outputs at the stage root.
    """
    logging.info(f"[correct-qc] run_qc regionalized={is_regionalized} n_regions={len(bed_files)}")
    qc_metrics: dict[str, dict] = {}

    for bed_file, region_tag in bed_files:
        # corrected bed written by flair correct
        corrected_bed = stage_dir / f"{region_tag}_all_corrected.bed"

        # Choose QC output directory
        region_qc_dir = (stage_dir / region_tag) if is_regionalized else stage_dir
        region_qc_dir.mkdir(parents=True, exist_ok=True)

        # Skip if corrected bed missing/empty (still emit a minimal TSV)
        if _bed_is_empty(corrected_bed):
            metrics = {"skipped_empty_bed": True, "n_corrected_reads": 0, "qc_runtime_sec": 0.0}
            write_metrics(region_qc_dir, "correct", metrics)
            qc_metrics[region_tag] = metrics
            logging.info(f"[correct-qc] Skipped: missing/empty corrected BED for {region_tag}: {corrected_bed}")
            continue

        try:
            metrics = collect(
                bed=corrected_bed,
                out_dir=region_qc_dir,
                n_input_reads=(None if is_regionalized else n_input_reads),
                align_sig=(None if is_regionalized else align_sig),
                genome_fa=genome_fa,
                runtime_sec=runtime_sec,
                regionalized=is_regionalized,
            )
            qc_metrics[region_tag] = metrics
        except Exception as e:
            logging.error(f"[correct-qc] QC failed for {region_tag}: {e}")

    return qc_metrics
