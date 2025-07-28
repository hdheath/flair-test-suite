# src/flair_test_suite/qc/correct_qc.py
# ----------------------------------
# QC collector for the 'correct' stage (revamped).
#  • Reports reads removed (count & %)
#  • Counts unique splice junctions before vs after correction
#  • Records both correct‑stage and QC runtimes

from __future__ import annotations
from pathlib import Path
import time

from . import register, write_metrics
from .qc_utils import count_lines, percent, count_unique_junctions, count_splice_junction_motifs
from ..lib.signature import qc_sidecar_path, load_marker
import json

__all__ = ["collect"]


@register("correct")
def collect(
    bed: Path,
    out_dir: Path,
    n_input_reads: int,
    align_sig: str,
    genome_fa: str,  # Path to genome FASTA for motif counting
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
    print(f"[DEBUG] correct_qc.collect genome_fa: {genome_fa}")

    # reads retained / removed
    n_corrected = count_lines(bed)
    n_removed = n_input_reads - n_corrected
    removed_pct = percent(n_removed, n_input_reads)

    # unique junctions after correction
    uniq_after = count_unique_junctions(bed)

    # unique junctions before correction (from align QC TSV)
    uniq_before = None
    try:
        run_dir = out_dir.parent.parent  # outputs/<run_id>/
        align_stage_dir = run_dir / "align" / align_sig
        qc_tsv = qc_sidecar_path(align_stage_dir, "align")
        print("QC TSV path:", qc_tsv)  # <-- Place this here
        if align_stage_dir.is_dir() and qc_tsv.exists():
            with open(qc_tsv) as fh:
                next(fh)  # skip header
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) == 2 and parts[0] == "unique_junctions":
                        uniq_before = int(parts[1])
                        break
    except Exception as e:
        print("Error reading align QC:", e)
        pass

    # Load motif counts before correction
    motif_counts_before = {}
    motif_json = align_stage_dir / "splice_site_motifs.json"
    try:
        if motif_json.exists():
            with open(motif_json) as fh:
                motif_counts_before = json.load(fh)
    except Exception as e:
        print("Error reading align motif counts JSON:", e)
        pass

    # Count motifs after correction
    try:
        motif_counts_after = count_splice_junction_motifs(
            bed_path=bed,
            fasta_path=Path(genome_fa),
            max_workers=4
        )
    except Exception as e:
        motif_counts_after = {}
        print(f"Warning: Splice motif counting (after) failed: {e}")

    # Compute diff and convert keys to strings
    motif_counts_after_str = {f"{k[0]}:{k[1]}" if isinstance(k, tuple) else str(k): v for k, v in motif_counts_after.items()}
    motif_counts_before_str = {f"{k[0]}:{k[1]}" if isinstance(k, tuple) else str(k): v for k, v in motif_counts_before.items()}
    motif_diff = {}
    for k in set(motif_counts_after_str) | set(motif_counts_before_str):
        motif_diff[k] = motif_counts_after_str.get(k, 0) - motif_counts_before_str.get(k, 0)
    motif_diff_str = {str(k): v for k, v in motif_diff.items()}

    # Save both before and after motif counts (and diff) to JSON
    motif_json_out = Path(out_dir) / "correct_splice_site_motifs.json"
    motif_json_data = {
        "after": motif_counts_after_str,
        "before": motif_counts_before_str,
        "diff": motif_diff_str,
    }
    with open(motif_json_out, "w") as fh:
        json.dump(motif_json_data, fh, indent=2)
    print(f"[DEBUG] Saved motifs before/after/diff to {motif_json_out}")

    # assemble QC metrics (motif counts NOT included)
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

