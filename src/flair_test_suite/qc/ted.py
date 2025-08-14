from __future__ import annotations

import io
import logging
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from . import register, write_metrics  # QC plumbing
from .qc_utils import count_lines

# If your CLI sets logging, you can remove/relax this.
logging.basicConfig(level=logging.DEBUG)

# ────────────────────────── tiny utils ──────────────────────────
def _which(exe: str) -> bool:
    from shutil import which
    return which(exe) is not None


def _resolve(p: Optional[str], data_dir: Path) -> Optional[Path]:
    if not p:
        return None
    pp = Path(p)
    return (data_dir / p) if not pp.is_absolute() else pp


def _run(cmd: List[str]) -> subprocess.CompletedProcess:
    logging.debug(f"[TED] RUN: {' '.join(cmd)}")
    try:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"[TED] Command failed ({' '.join(cmd)}): {e.stderr.strip()}")
        raise


def _read_bed6(path: Path) -> pd.DataFrame:
    """
    Read peaks as BED6 if possible; tolerate BED3 by adding Strand='.' (strand-insensitive).
    """
    try:
        return pd.read_csv(
            path, sep="\t", header=None, usecols=[0, 1, 2, 5],
            names=["Chrom", "Start", "End", "Strand"],
            dtype={"Chrom": str, "Start": np.int32, "End": np.int32, "Strand": str},
            comment="#",
        )
    except Exception:
        # Fallback: BED3
        logging.warning(f"[TED] {path} does not look like BED6; treating as BED3 (strand='.') which may zero-out recall.")
        df = pd.read_csv(
            path, sep="\t", header=None, usecols=[0, 1, 2],
            names=["Chrom", "Start", "End"],
            dtype={"Chrom": str, "Start": np.int32, "End": np.int32},
            comment="#",
        )
        df["Strand"] = "."
        return df


def _prepare_bed6_sorted(path: Path, tmps: List[Path]) -> Path:
    trimmed = path.with_suffix(path.suffix + ".trimmed.tmp")
    sorted_p = path.with_suffix(path.suffix + ".sorted.tmp")

    # Trim to first six columns of BED and write to temporary file
    res = _run(["cut", "-f", "1-6", str(path)])
    trimmed.write_text(res.stdout)

    # Sort the trimmed BED and capture output
    res = _run(["bedtools", "sort", "-i", str(trimmed)])
    sorted_p.write_text(res.stdout)

    try:
        trimmed.unlink()
    except Exception:
        pass
    tmps.append(sorted_p)
    logging.debug(f"[TED] Prepared sorted BED6: {sorted_p}")
    return sorted_p


def _run_closest(a: Path, b: Path) -> pd.DataFrame:
    res = _run(["bedtools", "closest", "-a", str(a), "-b", str(b), "-s", "-d"])
    return pd.read_csv(io.StringIO(res.stdout), sep="\t", header=None, comment="#")


def _extract_distance_and_peak(df: pd.DataFrame, label: str, max_dist: int) -> pd.DataFrame:
    orig_cols = df.shape[1]
    tx_series = df.iloc[:, 3].astype(str)
    dist_series = pd.to_numeric(df.iloc[:, orig_cols - 1], errors="coerce").fillna(max_dist + 1).astype(int)
    # Peak tuple (Chrom, Start, End) for -b
    peak_series = (df.iloc[:, 6].astype(str) + "_" + df.iloc[:, 7].astype(str) + "_" + df.iloc[:, 8].astype(str))
    out = pd.DataFrame({"isoform_id": tx_series, f"{label}_dist": dist_series, f"{label}_peak": peak_series})
    return out.groupby("isoform_id", as_index=False).first()


def _vectorized_overlap_counts(bed_df: pd.DataFrame, peaks_df: pd.DataFrame, window: int) -> Tuple[int, int]:
    """Return (#peaks matched by ≥1 isoform within window, total peaks)."""
    consumed = np.zeros(len(peaks_df), dtype=bool)
    for (c, s), chunk in bed_df.groupby(["Chrom", "Strand"], sort=False):
        mask = (peaks_df["Chrom"] == c) & (peaks_df["Strand"] == s)
        idxs = np.flatnonzero(mask)
        if idxs.size == 0:
            continue
        sub = peaks_df.iloc[idxs]
        starts = sub["Start"].to_numpy()
        ends = sub["End"].to_numpy()
        order = np.argsort(starts)
        starts = starts[order]
        ends = ends[order]
        global_i = idxs[order]
        for _, row in chunk.iterrows():
            lo, hi = int(row["Start"]) - window, int(row["End"]) + window
            r = np.searchsorted(starts, hi, side="right")
            if r == 0:
                continue
            hits = np.nonzero(ends[:r] >= lo)[0]
            if hits.size == 0:
                continue
            matched = global_i[hits]
            new_hits = matched[~consumed[matched]]
            consumed[new_hits] = True
    return int(consumed.sum()), int(len(peaks_df))


def _safe_f1(p: Optional[float], r: Optional[float]) -> Optional[float]:
    if p is None or r is None:
        return None
    s = p + r
    return (2 * p * r / s) if s > 0 else 0.0


def _extract_gene_id(name: str) -> Optional[str]:
    """
    Isoform name like: <read>_ENSG00000310526.1 → return ENSG00000310526 (drop version).
    If no recognizable gene suffix, return None.
    """
    if not name:
        return None
    tail = name.split("_")[-1]
    if tail.startswith(("ENS", "ENSG", "ENST", "ENSMUSG", "ENSMUST")):
        return tail.split(".")[0]
    m = re.search(r"(ENS[GMT][A-Z]*\d+)(?:\.\d+)?", name)
    return m.group(1) if m else None


def _isoform_counts(iso_bed: Path) -> Tuple[int, int]:
    df = pd.read_csv(iso_bed, sep="\t", header=None, comment="#", usecols=[3], names=["Name"], dtype=str)
    n_iso = len(df)
    genes = {g for g in (_extract_gene_id(nm) for nm in df["Name"].tolist()) if g}
    return n_iso, len(genes)


def _read_map_unique_reads(map_path: Path) -> int:
    """Count unique read IDs across all isoforms in *.isoform.read.map.txt."""
    if not map_path.exists() or map_path.stat().st_size == 0:
        return 0
    uniq = set()
    with open(map_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            _, rhs = line.split("\t", 1)
            for rid in rhs.split(","):
                rid = rid.strip()
                if rid:
                    uniq.add(rid)
    return len(uniq)


def _audit(
    rows: List[Dict],
    context: str,
    tag: str,
    label: str,
    path: Optional[Path],
    *,
    line_counter=count_lines,
    lines: Optional[int] = None,
    derived: Optional[Dict] = None,
    note: str = "",
) -> None:
    """Append a standardized audit record for a file."""
    exists = bool(path and path.exists())
    size = path.stat().st_size if exists else None
    if exists:
        if lines is not None:
            lines_est = lines
        elif line_counter:
            try:
                lines_est = line_counter(path)
            except Exception:
                lines_est = None
        else:
            lines_est = None
    else:
        lines_est = None
    rows.append(
        {
            "context": context,
            "tag": tag,
            "label": label,
            "path": str(path) if path else None,
            "exists": exists,
            "size_bytes": size,
            "lines_est": lines_est,
            "derived_count": derived,
            "note": note,
        }
    )


# ────────────────────────── discovery helpers ──────────────────────────
def _dirs_with_file(root: Path, stage: str, fname: str) -> List[Path]:
    base = root / stage
    if not base.exists():
        return []
    dirs = [p for p in base.iterdir() if p.is_dir() and (p / fname).exists()]
    dirs.sort(key=lambda p: (p / fname).stat().st_mtime, reverse=True)
    logging.debug(f"[TED] Searching {base} for {fname} -> {[str(d) for d in dirs]}")
    return dirs


def _build_region_metrics_index(run_root: Path) -> Dict[str, Dict]:
    """
    Build index: region_tag -> {'gene_count', 'transcript_count', 'reg_dir'}
    by scanning ALL regionalize/*/region_metrics.tsv.
    """
    idx: Dict[str, Dict] = {}
    regionalize_root = run_root / "regionalize"
    if not regionalize_root.exists():
        return idx
    for d in regionalize_root.iterdir():
        if not d.is_dir():
            continue
        metrics = d / "region_metrics.tsv"
        if not metrics.exists():
            continue
        try:
            df = pd.read_csv(metrics, sep="\t")
        except Exception as e:
            logging.warning(f"[TED] Failed reading {metrics}: {e}")
            continue
        if "region_tag" not in df.columns:
            continue
        for _, row in df.iterrows():
            tag = str(row["region_tag"])
            rec = {
                "gene_count": int(row.get("gene_count", 0)) if pd.notna(row.get("gene_count", np.nan)) else 0,
                "transcript_count": int(row.get("transcript_count", 0)) if pd.notna(row.get("transcript_count", np.nan)) else 0,
                "reg_dir": d,
            }
            idx[tag] = rec
    logging.debug(f"[TED] Region index entries: {len(idx)}")
    return idx


def _find_correct_bed_for_tag(run_root: Path, tag: str) -> Optional[Path]:
    for d in _dirs_with_file(run_root, "correct", f"{tag}_all_corrected.bed"):
        return d / f"{tag}_all_corrected.bed"
    return None


def _find_correct_bed_single(run_root: Path, run_id: str) -> Optional[Path]:
    for d in _dirs_with_file(run_root, "correct", f"{run_id}_all_corrected.bed"):
        return d / f"{run_id}_all_corrected.bed"
    return None


def _find_align_bam(run_root: Path, run_id: str) -> Optional[Path]:
    for d in _dirs_with_file(run_root, "align", f"{run_id}_flair.bam"):
        return d / f"{run_id}_flair.bam"
    return None

def _count_primary_alignments_bam(bam: Path) -> int:
    if not _which("samtools"):
        raise RuntimeError("TED requires 'samtools' on PATH to count primary alignments from BAM.")
    if not bam.exists() or bam.stat().st_size == 0:
        return 0
    # exclude unmapped(0x4), secondary(0x100), supplementary(0x800)
    res = _run(["samtools", "view", "-c", "-F", "2308", str(bam)])
    try:
        return int(res.stdout.strip() or "0")
    except Exception:
        return 0


def _is_region_tag(stem: str) -> bool:
    """chr_start_end[...].isoforms.bed → True if start/end are integers."""
    m = re.match(r"^(.+?)_(\d+)_(\d+)$", stem)  # stem has no .isoforms suffix when we call this
    return m is not None


def _tss_tts_metrics_full(iso_bed: Path, peaks: Dict[str, Optional[Path]], window: int,
                          audit_rows: List[Dict], tag_ctx: str) -> Dict[str, Optional[float]]:
    """
    Compute precision/recall/F1 for experimental and reference TSS/TTS:
      - experimental: keys 'prime5' → 5prime_*, 'prime3' → 3prime_*
      - reference:    keys 'ref_prime5' → ref5prime_*, 'ref_prime3' → ref3prime_*
    If a file is missing, that side returns None for all three metrics.
    """
    have_any = any(peaks.values())
    if have_any and not _which("bedtools"):
        raise RuntimeError("TED requires 'bedtools' on PATH when TSS/TTS files are provided.")

    full_bed_df = pd.read_csv(
        iso_bed, sep="\t", header=None, comment="#",
        usecols=[0, 1, 2, 3, 5],
        names=["Chrom", "Start", "End", "Name", "Strand"],
        dtype={"Chrom": str, "Start": np.int32, "End": np.int32, "Name": str, "Strand": str},
    )
    n_tx = len(full_bed_df)
    bed_df = full_bed_df[["Chrom", "Start", "End", "Strand"]]

    tmp_local: List[Path] = []
    metrics: Dict[str, Optional[float]] = {
        "5prime_precision": None, "5prime_recall": None, "5prime_f1": None,
        "3prime_precision": None, "3prime_recall": None, "3prime_f1": None,
        "ref5prime_precision": None, "ref5prime_recall": None, "ref5prime_f1": None,
        "ref3prime_precision": None, "ref3prime_recall": None, "ref3prime_f1": None,
    }
    if not have_any:
        return metrics

    try:
        bed_sorted = _prepare_bed6_sorted(iso_bed, tmp_local)

        def _side(key_cfg: str, label: str) -> Tuple[Optional[float], Optional[float], Optional[float]]:
            pth = peaks.get(key_cfg)
            if pth is None or not pth.exists() or pth.stat().st_size == 0:
                note = "missing" if (pth is None or not pth.exists()) else "empty"
                _audit(audit_rows, tag_ctx, label, f"{key_cfg}_peaks", pth, note=note)
                return None, None, None

            # Sort/trim peaks, too
            peaks_sorted = _prepare_bed6_sorted(pth, tmp_local)
            dfc = _run_closest(bed_sorted, peaks_sorted)
            ddf = _extract_distance_and_peak(dfc, label, window)
            m = int((ddf[f"{label}_dist"] <= window).sum())
            precision = (m / n_tx) if n_tx else None

            peaks_df = _read_bed6(pth)
            c, t = _vectorized_overlap_counts(bed_df, peaks_df, window)
            recall = (c / t) if t else None
            f1 = _safe_f1(precision, recall)

            _audit(
                audit_rows,
                tag_ctx,
                label,
                f"{key_cfg}_peaks",
                pth,
                derived={"matches_le_window": m, "peaks_total": t},
                note=f"precision={precision}, recall={recall}, f1={f1}",
            )
            return precision, recall, f1

        p5, r5, f5 = _side("prime5", "5")
        p3, r3, f3 = _side("prime3", "3")
        rp5, rr5, rf5 = _side("ref_prime5", "ref5")
        rp3, rr3, rf3 = _side("ref_prime3", "ref3")

        metrics.update({
            "5prime_precision": p5, "5prime_recall": r5, "5prime_f1": f5,
            "3prime_precision": p3, "3prime_recall": r3, "3prime_f1": f3,
            "ref5prime_precision": rp5, "ref5prime_recall": rr5, "ref5prime_f1": rf5,
            "ref3prime_precision": rp3, "ref3prime_recall": rr3, "ref3prime_f1": rf3,
        })
        return metrics
    finally:
        for f in tmp_local:
            try:
                f.unlink()
            except Exception:
                pass


# ────────────────────────── main collector ──────────────────────────
@register("TED")
@register("collapse")
def collect(stage_dir: Path, cfg) -> None:
    """
    Writes:
      - TED.tsv: per-region rows if regionalized else one row
      - TED.audit.tsv: everything we touched + counts
      - TED sidecar with simple totals (no means)
    Implements:
      * region_genes_expected / region_isoforms_expected via an index over ALL regionalize dirs
      * assigned_pct denominator: 
          - collapse: corrected BED line count (per-region if regionalized, else single)
          - transcriptome: primary alignments in BAM (per-region regional BAM, else full align BAM)
    """
    logging.debug("[TED] Starting TED QC collection.")
    stage_name = stage_dir.parent.name           # collapse | transcriptome
    run_root = stage_dir.parent.parent           # .../outputs/<run_id>
    run_id = run_root.name
    logging.debug(f"[TED] Stage name: {stage_name}, Run ID: {run_id}")

    # cfg access helper (object or dict)
    def _cfg_get(path: List[str], default=None):
        cur = cfg
        try:
            for k in path:
                cur = getattr(cur, k)
            return cur
        except Exception:
            cur = cfg
            for k in path:
                cur = cur.get(k, {}) if isinstance(cur, dict) else {}
            return cur if cur != {} else default

    data_dir = Path(_cfg_get(["run", "data_dir"], "."))
    logging.debug(f"[TED] Data directory set to: {data_dir}")

    # Optional QC block
    qc_block = _cfg_get(["qc", stage_name, "TED"], {}) or {}
    window = int(qc_block.get("window", 50))
    logging.debug(f"[TED] Window: {window}")
    peaks_cfg = {
        "prime5": qc_block.get("experiment_5_prime_regions_bed_file"),
        "prime3": qc_block.get("experiment_3_prime_regions_bed_file"),
        "ref_prime5": qc_block.get("reference_5_prime_regions_bed_file"),
        "ref_prime3": qc_block.get("reference_3_prime_regions_bed_file"),
    }
    logging.debug(f"[TED] Peaks configuration (raw): {peaks_cfg}")

    # Build region metrics index across ALL regionalize runs
    reg_index = _build_region_metrics_index(run_root)

    rows: List[Dict] = []
    audit_rows: List[Dict] = []

    # Determine if this is regionalized by matching chr_start_end pattern
    iso_beds = sorted(stage_dir.glob("*.isoforms.bed"))
    regional_files = [p for p in iso_beds if _is_region_tag(p.stem.replace(".isoforms", ""))]
    is_regionalized = len(regional_files) > 0
    logging.debug(f"[TED] Isoform beds: {[p.name for p in iso_beds]}")
    logging.debug(f"[TED] Regional isoform beds: {[p.name for p in regional_files]}")
    logging.debug(f"[TED] is_regionalized: {is_regionalized}")

    if is_regionalized:
        # ── Regionalized: iterate over discovered isoform files
        for iso_bed in regional_files:
            tag = iso_bed.stem.replace(".isoforms", "")
            logging.debug(f"[TED] Processing region tag: {tag}")
            try:
                chrom, start, end = tag.split("_", 3)
                start_i, end_i = int(start), int(end)
            except Exception:
                logging.warning(f"[TED] Could not parse region tag from {iso_bed.name}; skipping")
                continue

            _audit(audit_rows, "region", tag, "isoforms_bed", iso_bed)

            map_txt = stage_dir / f"{tag}.isoform.read.map.txt"
            if map_txt.exists():
                _audit(audit_rows, "region", tag, "map_txt", map_txt)
            else:
                _audit(audit_rows, "region", tag, "map_txt", map_txt, note="missing")

            n_iso, n_genes_obs = _isoform_counts(iso_bed)
            assigned_reads = _read_map_unique_reads(map_txt)
            logging.debug(f"[TED] Region {tag} - n_iso: {n_iso}, n_genes_obs: {n_genes_obs}, assigned_reads: {assigned_reads}")

            # expected region metrics via index
            reg_rec = reg_index.get(tag, {})
            reg_gene_cnt = int(reg_rec.get("gene_count", 0) or 0)
            reg_tx_cnt   = int(reg_rec.get("transcript_count", 0) or 0)
            reg_dir_for_tag: Optional[Path] = reg_rec.get("reg_dir")
            logging.debug(f"[TED] For region {tag}, expected genes: {reg_gene_cnt}, transcripts: {reg_tx_cnt}")

            # peaks: prefer sliced if available
            peaks = {}
            for key, conf in peaks_cfg.items():
                if not conf:
                    peaks[key] = None
                    logging.debug(f"[TED] No configuration for peaks key '{key}'")
                    continue
                base = Path(conf).name
                sliced = (reg_dir_for_tag / f"{tag}_{base}") if reg_dir_for_tag else None
                if sliced and sliced.exists() and sliced.stat().st_size > 0:
                    peaks[key] = sliced
                else:
                    rp = _resolve(conf, data_dir)
                    peaks[key] = rp if (rp and rp.exists() and rp.stat().st_size > 0) else None

            logging.debug(f"[TED] For region {tag}, resolved peaks: { {k: str(v) if v else None for k,v in peaks.items()} }")

            tss_tts = _tss_tts_metrics_full(iso_bed, peaks, window, audit_rows, f"region:{tag}")
            logging.debug(f"[TED] Region {tag} TSS/TTS metrics: {tss_tts}")

            # Denominator per rule
            if stage_name == "collapse":
                corr_bed = _find_correct_bed_for_tag(run_root, tag)
                if corr_bed is None:
                    logging.warning(f"[TED] Corrected BED for {tag} not found; assigned_pct will be None")
                    denom = None
                    _audit(audit_rows, "region", tag, "corrected_bed", None, note="missing")
                else:
                    denom = count_lines(corr_bed)
                    _audit(
                        audit_rows,
                        "region",
                        tag,
                        "corrected_bed",
                        corr_bed,
                        lines=denom,
                        derived={"denom": denom},
                    )
                    logging.debug(f"[TED] For region {tag}, corrected BED: {corr_bed.name}, lines: {denom}")
            else:
                reg_bam = (reg_dir_for_tag / f"{tag}.bam") if reg_dir_for_tag else None
                if reg_bam and reg_bam.exists():
                    denom = _count_primary_alignments_bam(reg_bam)
                    _audit(
                        audit_rows,
                        "region",
                        tag,
                        "regional_bam",
                        reg_bam,
                        line_counter=None,
                        derived={"primary_alignments": denom},
                    )
                else:
                    logging.warning(f"[TED] Regional BAM for {tag} not found; assigned_pct will be None")
                    denom = None
                    _audit(
                        audit_rows,
                        "region",
                        tag,
                        "regional_bam",
                        reg_bam,
                        line_counter=None,
                        note="missing",
                    )

            row = {
                "run_id": run_id,
                "stage": stage_name,
                "region_tag": tag,
                "chrom": chrom, "start": start_i, "end": end_i,
                "span_bp": int(end_i - start_i + 1),
                "region_genes_expected": reg_gene_cnt,
                "region_isoforms_expected": reg_tx_cnt,
                "isoforms_observed": n_iso,
                "genes_observed": n_genes_obs,
                "assigned_primary_reads": int(assigned_reads),
                "align_primary_total": int(denom) if (denom is not None) else None,
                "assigned_pct": (assigned_reads / denom) if (denom and denom > 0) else None,
                **tss_tts,
            }
            logging.debug(f"[TED] Region {tag} final row: {row}")
            if reg_dir_for_tag is None:
                logging.debug(f"[TED] No regionalize index entry for {tag}; expected metrics defaulted to 0")
            rows.append(row)
    else:
        # ── Single (non-regionalized): one row ──
        iso_bed = stage_dir / f"{run_id}.isoforms.bed"
        if not iso_bed.exists() or iso_bed.stat().st_size == 0:
            raise RuntimeError(f"[TED] Isoforms BED not found or empty: {iso_bed}")
        map_txt = stage_dir / f"{run_id}.isoform.read.map.txt"

        _audit(audit_rows, "single", run_id, "isoforms_bed", iso_bed)

        if map_txt.exists():
            _audit(audit_rows, "single", run_id, "map_txt", map_txt)
        else:
            _audit(audit_rows, "single", run_id, "map_txt", map_txt, note="missing")

        n_iso, n_genes_obs = _isoform_counts(iso_bed)
        assigned_reads = _read_map_unique_reads(map_txt)
        logging.debug(f"[TED] Single mode - n_iso: {n_iso}, n_genes_obs: {n_genes_obs}, assigned_reads: {assigned_reads}")

        # peaks: global only
        peaks = {}
        for key, conf in peaks_cfg.items():
            if not conf:
                peaks[key] = None
                continue
            rp = _resolve(conf, data_dir)
            if rp and rp.exists() and rp.stat().st_size > 0:
                peaks[key] = rp
            else:
                peaks[key] = None
        logging.debug(f"[TED] Single mode peaks resolved: { {k: str(v) if v else None for k, v in peaks.items()} }")
        tss_tts = _tss_tts_metrics_full(iso_bed, peaks, window, audit_rows, "single")

        if stage_name == "collapse":
            corr_bed = _find_correct_bed_single(run_root, run_id)
            if corr_bed:
                denom = count_lines(corr_bed)
                _audit(
                    audit_rows,
                    "single",
                    run_id,
                    "corrected_bed",
                    corr_bed,
                    lines=denom,
                    derived={"denom": denom},
                )
                logging.debug(f"[TED] Found single corrected BED: {corr_bed.name} with {denom} lines")
            else:
                logging.warning("[TED] Single corrected BED not found; assigned_pct will be None")
                denom = None
                _audit(audit_rows, "single", run_id, "corrected_bed", None, note="missing")
        else:
            aln_bam = _find_align_bam(run_root, run_id)
            if aln_bam:
                denom = _count_primary_alignments_bam(aln_bam)
                _audit(
                    audit_rows,
                    "single",
                    run_id,
                    "align_bam",
                    aln_bam,
                    line_counter=None,
                    derived={"primary_alignments": denom},
                )
                logging.debug(f"[TED] Found single align BAM: {aln_bam.name}")
            else:
                logging.warning("[TED] Align BAM not found; assigned_pct will be None")
                denom = None
                _audit(
                    audit_rows,
                    "single",
                    run_id,
                    "align_bam",
                    None,
                    line_counter=None,
                    note="missing",
                )

        row = {
            "run_id": run_id,
            "stage": stage_name,
            "region_tag": None,
            "chrom": None, "start": None, "end": None, "span_bp": None,
            "region_genes_expected": None,
            "region_isoforms_expected": None,
            "isoforms_observed": n_iso,
            "genes_observed": n_genes_obs,
            "assigned_primary_reads": int(assigned_reads),
            "align_primary_total": int(denom) if (denom is not None) else None,
            "assigned_pct": (assigned_reads / denom) if (denom and denom > 0) else None,
            **tss_tts,
        }
        logging.debug(f"[TED] Single mode final row: {row}")
        rows.append(row)

    # write TSV
    out_tsv = stage_dir / "TED.tsv"
    df = pd.DataFrame(rows)
    if is_regionalized:
        df.to_csv(out_tsv, sep="\t", index=False)
    else:
        drop_cols = [
            "region_tag",
            "chrom",
            "start",
            "end",
            "span_bp",
            "region_genes_expected",
            "region_isoforms_expected",
        ]
        df.drop(columns=drop_cols, inplace=True)
        df.to_csv(out_tsv, sep="\t", index=False)
    logging.info(f"[TED] Wrote TED.tsv with {len(rows)} rows at {out_tsv}")

    # write AUDIT TSV
    audit_tsv = stage_dir / "TED.audit.tsv"
    pd.DataFrame(audit_rows, columns=[
        "context", "tag", "label", "path", "exists", "size_bytes", "lines_est", "derived_count", "note"
    ]).to_csv(audit_tsv, sep="\t", index=False)
    logging.info(f"[TED] Wrote audit log at {audit_tsv}")

    # sidecar summary (simple totals—no means)
    assigned_total = int(sum((r.get("assigned_primary_reads") or 0) for r in rows))
    metrics = {
        "window": int(window),
        "rows": len(rows),
        "assigned_primary_reads_total": assigned_total,
        "tsv": str(out_tsv),
        "audit": str(audit_tsv),
    }
    write_metrics(stage_dir, "TED", metrics)
    logging.info(f"[TED] Sidecar metrics written: {metrics}")
