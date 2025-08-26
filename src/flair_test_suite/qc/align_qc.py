# src/flair_test_suite/qc/align_qc.py
# ---------------------------------
# Enhanced QC collector for the 'align' stage.
# Requirements:
#  • pysam            # for parsing BAM alignments
#  • samtools         # for generating alignment statistics
#  • matplotlib       # for plotting histograms (bundled in flair-test-suite env)
#  • tempfile, subprocess, time, statistics, json
# Summary:
#   This module gathers and reports multiple quality-control metrics after FLAIR's
#   `align` stage. Key outputs include:
#     • TSV metrics file (align_qc.tsv)
#     • Three PNG histograms (MAPQ, read identity, read length)
#     • JSON file of splice junction motif counts (splice_site_motifs.json)
#   Metrics computed:
#     - Number and percentage of reads retained (BED lines)
#     - Mapping quality distribution
#     - Read identity and soft-clipped read percentage
#     - Unique splice junction count
#     - Alignment and QC runtimes
# Functions:
#   collect(...)
#     Main function: executes all QC steps, writes metrics, plots, and JSON.
#   _hist(vals, filename, xlabel)
#     Internal helper: generates and saves a histogram for a list of values.

from __future__ import annotations
from pathlib import Path
import json
import subprocess
import tempfile
import time
from statistics import mean, median

# Import utilities for QC: sidecar path, marker loading, registration, metrics write
import pysam  # noqa: E402
from ..plotting.align_histograms import generate_histograms

from . import register, write_metrics
from .qc_utils import (
    count_lines,
    percent,
    iter_primary,
    count_unique_junctions,
    SAMPLE_LIMIT,
    count_splice_junction_motifs
)

__all__ = ["collect"]

@register("align")
def collect(
    bam: Path,
    out_dir: Path,
    n_input_reads: int | None,
    genome_fa: str,
    runtime_sec: float | None = None,
    read_count_method: str | None = None,
) -> dict:
    """
    Main QC collector for the 'align' stage.

    Arguments:
      bam           : Path to aligned BAM file
      out_dir       : Directory to save QC outputs
      n_input_reads : Optional number of input reads; when provided a mapped
                      percentage is calculated
      genome_fa     : Reference FASTA used for motif counting
      runtime_sec   : Time taken by the align stage (seconds)
      read_count_method : "exact" or "estimated" depending on counting strategy

    Returns:
      Dictionary of collected QC metrics.
    """
    qc_start = time.time()

    # 1. Count retained reads and compute mapped percentage
    bed = bam.with_suffix(".bed")
    retained = count_lines(bed)
    mapped_pct = percent(retained, n_input_reads) if n_input_reads else None

    # 2. Extract MAPQ and read length distributions via samtools stats
    mapq_vals: list[int] = []
    read_len_vals: list[int] = []
    with tempfile.TemporaryDirectory() as tmpd:
        stats_out = Path(tmpd) / "stats.txt"
        try:
            subprocess.run(
                ["samtools", "stats", "-F", "0x904", "-o", str(stats_out), str(bam)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
        except subprocess.CalledProcessError:
            with open(stats_out, "w") as fh:
                subprocess.run(
                    ["samtools", "stats", "-F", "0x904", str(bam)],
                    check=True, stdout=fh, stderr=subprocess.DEVNULL
                )
        for line in stats_out.read_text().splitlines():
            if line.startswith("MAPQ\t"):
                _, v, c = line.split("\t")
                mapq_vals.extend([int(v)] * int(c))
            elif line.startswith("RL\t"):
                _, l, c = line.split("\t")
                read_len_vals.extend([int(l)] * int(c))

    # 3. Sample alignments to compute identity and soft-clip stats
    identity_vals: list[float] = []
    softclip_n = 0
    total_sampled = 0
    with pysam.AlignmentFile(bam, "rb") as bam_f:
        for aln in iter_primary(bam_f, SAMPLE_LIMIT):
            total_sampled += 1
            if aln.has_tag("NM") and aln.query_length:
                identity_vals.append(1 - aln.get_tag("NM") / aln.query_length)
            if any(op == 4 for op, _ in (aln.cigartuples or [])):
                softclip_n += 1
    softclip_pct = percent(softclip_n, total_sampled)

    # 4. Count unique splice junctions
    unique_juncs = count_unique_junctions(bed)

    # 5. Plot histograms and save filenames under <stage>/qc
    qc_dir = Path(out_dir) / "qc"
    mapping = generate_histograms(
        mapq_vals,
        [round(v * 100, 1) for v in identity_vals],
        read_len_vals,
        qc_dir,
    )
    mapq_png = mapping.get("mapq")
    id_png = mapping.get("identity")
    len_png = mapping.get("length")

    # 6. Count splice junction motifs (4-mer) and write JSON
    try:
        motif_counts = count_splice_junction_motifs(
            bed_path=bed,
            fasta_path=Path(genome_fa),
            max_workers=4
        )
    except Exception as e:
        motif_counts = {}
        print(f"Warning: Splice junction motif counting failed: {e}")
    motif_counts_str = {f"{k[0]}:{k[1]}": v for k, v in motif_counts.items()}
    with open(qc_dir / "splice_site_motifs.json", "w") as fh:
        json.dump(motif_counts_str, fh, indent=2)

    # 7. Compile metrics and write outputs
    metrics = {
        "n_input_reads":      n_input_reads,
        "read_count_method":  read_count_method or "unknown",
        "n_retained_bed":     retained,
        "mapped_pct":         mapped_pct,
        "mean_identity":      round(mean(identity_vals)*100, 2) if identity_vals else 0.0,
        "median_identity":    round(median(identity_vals)*100, 2) if identity_vals else 0.0,
        "n_softclip":         softclip_n,
        "softclip_pct":       softclip_pct,
        "unique_junctions":   unique_juncs,
        "align_runtime_sec":  round(runtime_sec, 2) if runtime_sec else None,
        "qc_runtime_sec":     round(time.time() - qc_start, 2)
    }
    write_metrics(out_dir, "align", metrics)
    with open(qc_dir / "align_plot_manifest.json", "w") as fh:
        json.dump({"mapq": mapq_png, "identity": id_png, "length": len_png}, fh, indent=2)

    return metrics
