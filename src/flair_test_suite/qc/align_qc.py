# src/flair_test_suite/qc/align_qc.py
from __future__ import annotations
import subprocess
from pathlib import Path
from . import register, write_metrics         #  ← relative import

@register("align")
def collect(
    bam: Path,
    out_dir: Path,
    n_input_reads: int,
    runtime_sec: float | None = None,
) -> dict:

    bed       = bam.with_suffix(".bed")
    retained  = int(subprocess.check_output(["wc", "-l", str(bed)], text=True).split()[0])
    mapped_pct = round(100 * retained / n_input_reads, 2)

    metrics = {
        "n_total_reads":  n_input_reads,
        "n_retained_bed": retained,
        "mapped_pct":     mapped_pct,
        "retained_pct":   mapped_pct,   # same definition
    }
    if runtime_sec is not None:
        metrics["runtime_sec"] = round(runtime_sec, 2)

    write_metrics(out_dir / "align_qc.tsv", metrics)
    return metrics
