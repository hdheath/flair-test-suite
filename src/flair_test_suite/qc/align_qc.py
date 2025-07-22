# src/flair_test_suite/qc/align_qc.py
# ---------------------------------
# QC collector for the 'align' stage.
# Computes read retention metrics and writes them to a TSV side‑car.

from __future__ import annotations
import subprocess      # to run shell commands (e.g., wc)
from pathlib import Path

# import the registry decorator and helper to write TSVs
from . import register, write_metrics


@register("align")  # registers this function under QC_REGISTRY['align']
def collect(
    bam: Path,
    out_dir: Path,
    n_input_reads: int,
    runtime_sec: float | None = None,
) -> dict:
    """
    QC collector for the align stage.

    Parameters:
      - bam: Path to the BAM output of flair align
      - out_dir: Directory where side‑cars should be written
      - n_input_reads: Number of reads counted before alignment
      - runtime_sec: Optional runtime in seconds (from StageBase)

    Returns:
      - metrics dict containing counts and percentages
    """
    # Derive the BED path by changing .bam → .bed
    bed = bam.with_suffix(".bed")

    # Count number of lines in the BED (each line = one mapped read)
    retained = int(
        subprocess.check_output(
            ["wc", "-l", str(bed)], text=True
        ).split()[0]
    )

    # Compute percentage of reads that mapped
    mapped_pct = round(100 * retained / n_input_reads, 2)

    # Prepare the metrics dictionary:
    metrics = {
        # original input count (from build_cmd)
        "n_input_reads":   n_input_reads,
        # number of retained reads in BED file
        "n_retained_bed":  retained,
        # percent mapped
        "mapped_pct":      mapped_pct,
        # same as mapped_pct, kept for naming consistency
        "retained_pct":    mapped_pct,
    }

    # If runtime was provided, include it rounded to 2 decimals
    if runtime_sec is not None:
        metrics["runtime_sec"] = round(runtime_sec, 2)

    # Write metrics out to <out_dir>/align_qc.tsv via the helper
    write_metrics(out_dir, "align", metrics)

    # Return the same dict so StageBase can embed it in .metadata
    return metrics

