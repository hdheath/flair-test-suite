# src/flair_test_suite/qc/correct_qc.py
# ----------------------------------
# QC collector for the 'correct' stage (revamped).
#  • Reports reads removed (count & %)
#  • Counts unique splice junctions before vs after correction
#  • Records both correct‑stage and QC runtimes

from __future__ import annotations
from pathlib import Path
import time

from . import register, write_metrics, qc_sidecar_path
from .qc_utils import count_lines, percent, count_unique_junctions

__all__ = ["collect"]


@register("correct")
def collect(
    bed: Path,
    out_dir: Path,
    n_input_reads: int,
    runtime_sec: float | None = None,
) -> dict:
    """QC collector for the correct stage.

    Reports:
      • n_input_reads, n_corrected_reads
      • n_removed_reads, removed_pct
      • unique_junc_before (from align QC)
      • unique_junc_after (recount)
      • correct_runtime_sec (stage runtime)
      • qc_runtime_sec (time to compute these QC metrics)
    """
    qc_start = time.time()

    # reads retained / removed
    n_corrected = count_lines(bed)
    n_removed = n_input_reads - n_corrected
    removed_pct = percent(n_removed, n_input_reads)

    # unique junctions after correction
    uniq_after = count_unique_junctions(bed)

    # unique junctions before correction (from align QC TSV)
    uniq_before = None
    try:
        run_dir = out_dir.parent  # outputs/<run_id>/
        align_stage_dir = next((p for p in (run_dir / "align").iterdir() if p.is_dir()), None)
        if align_stage_dir:
            qc_tsv = qc_sidecar_path(align_stage_dir, "align")
            if qc_tsv.exists():
                with open(qc_tsv) as fh:
                    for line in fh:
                        if line.startswith("unique_junctions\t"):
                            uniq_before = int(line.split("\t")[1])
                            break
    except Exception:
        # leave uniq_before as None on any error
        pass

    # assemble QC metrics
    metrics = {
        "n_input_reads":      n_input_reads,
        "n_corrected_reads":  n_corrected,
        "n_removed_reads":    n_removed,
        "removed_pct":        removed_pct,
        "unique_junc_before": uniq_before,
        "unique_junc_after":  uniq_after,
        "correct_runtime_sec": round(runtime_sec, 2) if runtime_sec is not None else None,
        "qc_runtime_sec":     round(time.time() - qc_start, 2),
    }

    write_metrics(out_dir, "correct", metrics)
    return metrics


