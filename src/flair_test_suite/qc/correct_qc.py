# src/flair_test_suite/qc/correct_qc.py
# ----------------------------------
# QC collector for the 'correct' stage.
# Computes retention metrics after correction and writes them to a TSV side‑car.

from __future__ import annotations
import subprocess      # to run shell commands (e.g., wc)
from pathlib import Path

# import the registry decorator and helper to write TSVs
from . import register, write_metrics


@register("correct")  # registers this function under QC_REGISTRY['correct']
def collect(
    bed: Path,
    out_dir: Path,
    n_input_reads: int,
    runtime_sec: float | None = None,
) -> dict:
    """
    QC collector for the correct stage.

    Parameters:
      - bed: Path to the BED output of flair correct
      - out_dir: Directory where side‑cars should be written
      - n_input_reads: Number of reads input into correct stage
      - runtime_sec: Optional runtime in seconds (from StageBase)

    Returns:
      - metrics dict containing counts and retention percentage
    """
    # Count number of lines in the BED (each line = one corrected read)
    retained = int(
        subprocess.check_output(
            ["wc", "-l", str(bed)], text=True
        ).split()[0]
    )

    # Compute percentage of reads retained after correction
    retained_pct = round(100 * retained / n_input_reads, 2)

    # Prepare the metrics dictionary:
    metrics = {
        # original input read count passed from align stage
        "n_input_reads":     n_input_reads,
        # number of reads in the corrected BED file
        "n_corrected_reads": retained,
        # percent retained after correction
        "retained_pct":      retained_pct,
    }

    # If runtime was provided, include it rounded to 2 decimals
    if runtime_sec is not None:
        metrics["runtime_sec"] = round(runtime_sec, 2)

    # Write metrics out to <out_dir>/correct_qc.tsv via the helper
    write_metrics(out_dir, "correct", metrics)

    # Return metrics dict so StageBase can embed it in .metadata
    return metrics

