# src/flair_test_suite/qc/correct_qc.py
# ----------------------------------
# QC collector for the 'correct' stage (revamped).
# Requirements:
#  • pysam            # for parsing BAM alignments (motif counting)
#  • tempfile, subprocess, time, statistics, json
#  • qc_utils helpers: count_lines, percent, count_unique_junctions, count_splice_junction_motifs
# Summary:
#   This module evaluates quality-control metrics after FLAIR's `correct` stage. Key outputs include:
#     • TSV metrics file (correct_qc.tsv)
#     • JSON file of splice junction motif counts before/after/diff
#   Metrics computed:
#     - Number and percentage of reads removed
#     - Unique splice junction counts before and after correction
#     - Stage and QC runtimes
# Functions:
#   collect(...)
#     Main function: computes QC metrics, reads upstream align QC data, counts motifs, writes outputs.

from __future__ import annotations
from pathlib import Path
import time
import json

# Register collector and write helper
from . import register, write_metrics
# Utilities for QC calculations
from .qc_utils import (
    count_lines,
    percent,
    count_unique_junctions,
    count_splice_junction_motifs
)
# For locating previous QC sidecar files
from ..lib.signature import qc_sidecar_path, load_marker

__all__ = ["collect"]

@register("correct")
def collect(
    bed: Path,
    out_dir: Path,
    n_input_reads: int,
    align_sig: str,
    genome_fa: str,  # Reference FASTA for motif extraction
    runtime_sec: float | None = None,
) -> dict:
    """
    QC collector for the 'correct' stage.

    Arguments:
      bed           : Path to corrected reads BED file
      out_dir       : Directory to save QC outputs
      n_input_reads : Number of reads input to correction (from align stage)
      align_sig     : Signature of the align run (to locate upstream QC)
      genome_fa     : Reference FASTA for motif counting
      runtime_sec   : Time taken by the correct stage (seconds)

    Returns:
      Dictionary of collected QC metrics.
    """
    # Start QC timer
    qc_start = time.time()

    # 1. Count reads after correction and compute removal metrics
    n_corrected = count_lines(bed)
    n_removed = n_input_reads - n_corrected
    removed_pct = percent(n_removed, n_input_reads)

    # 2. Unique splice junctions after correction
    uniq_after = count_unique_junctions(bed)

    # 3. Load unique junctions from before correction (align QC TSV)
    uniq_before = None
    try:
        # Construct path to align QC TSV using run structure
        run_dir = out_dir.parent.parent  # outputs/<run_id>/<stage_sig>/correct -> outputs/<run_id>
        align_stage_dir = run_dir / "align" / align_sig
        qc_tsv = qc_sidecar_path(align_stage_dir, "align")
        if align_stage_dir.is_dir() and qc_tsv.exists():
            with open(qc_tsv) as fh:
                next(fh)  # skip header line
                for line in fh:
                    key, val = line.rstrip("\n").split("\t")
                    if key == "unique_junctions":
                        uniq_before = int(val)
                        break
    except Exception:
        # If reading fails, leave uniq_before as None
        pass

    # 4. Load splice site motif counts before correction (JSON)
    motif_counts_before: dict = {}
    try:
        motif_json_path = align_stage_dir / "splice_site_motifs.json"
        if motif_json_path.exists():
            with open(motif_json_path) as fh:
                motif_counts_before = json.load(fh)
    except Exception:
        pass

    # 5. Count splice site motifs after correction
    try:
        raw_after = count_splice_junction_motifs(
            bed_path=bed,
            fasta_path=Path(genome_fa),
            max_workers=4
        )
        # Convert tuple keys to string "chr:pos"
        motif_counts_after = {f"{k[0]}:{k[1]}": v for k, v in raw_after.items()}
    except Exception:
        motif_counts_after = {}

    # 6. Compute motif diff (after - before)
    motif_diff: dict = {}
    for k in set(motif_counts_after) | set(motif_counts_before):
        motif_diff[k] = motif_counts_after.get(k, 0) - motif_counts_before.get(k, 0)

    # 7. Save motif JSON with before/after/diff
    out_json = Path(out_dir) / "correct_splice_site_motifs.json"
    with open(out_json, "w") as fh:
        json.dump({
            "before": motif_counts_before,
            "after": motif_counts_after,
            "diff": motif_diff,
        }, fh, indent=2)

    # 8. Assemble QC metrics
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

    # Write metrics TSV
    write_metrics(out_dir, "correct", metrics)
    return metrics

