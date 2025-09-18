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
import logging
import os

logger = logging.getLogger(__name__)

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

    # Read tunables from environment (with sensible defaults)
    stats_threads = int(os.getenv("FTS_QC_STATS_THREADS", "4"))
    motif_limit = int(os.getenv("FTS_QC_MAX_BED_LINES", "1"))
    junc_limit = int(os.getenv("FTS_QC_JUNC_LIMIT", "0"))  # 0 => no cap
    sample_limit = int(os.getenv("FTS_QC_SAMPLE_LIMIT", str(SAMPLE_LIMIT)))
    motif_workers = int(os.getenv("FTS_QC_CPUS", "4"))
    junc_mode = os.getenv("FTS_QC_JUNC_MODE", "all").strip().lower()  # all|sample|skip
    stats_mode = os.getenv("FTS_QC_STATS_MODE", "all").strip().lower()  # all|sample|skip
    stats_timeout_env = os.getenv("FTS_QC_STATS_TIMEOUT", "").strip()
    try:
        stats_timeout = int(stats_timeout_env) if stats_timeout_env else None
    except Exception:
        stats_timeout = None
    logger.info(
        f"[align_qc] Tunables: stats_threads={stats_threads}, sample_limit={sample_limit}, "
        f"junc_limit={junc_limit or 'none'}, motif_limit={motif_limit}, motif_workers={motif_workers}, "
        f"stats_mode={stats_mode}, stats_timeout={stats_timeout or 'none'}"
    )

    # 1. Count retained reads and compute mapped percentage
    bed = bam.with_suffix(".bed")
    logger.info(f"[align_qc] Inputs: BAM={bam}, BED={bed}, genome={genome_fa}")
    retained = count_lines(bed)
    mapped_pct = percent(retained, n_input_reads) if n_input_reads else None
    logger.info(f"[align_qc] Retained BED lines: {retained}; n_input_reads={n_input_reads}")

    # 2. MAPQ and read length distributions: choose 'all' (samtools), 'sample' (from iter), or 'skip'
    mapq_vals: list[int] = []
    read_len_vals: list[int] = []
    need_sample_for_stats = (stats_mode != "all")
    if stats_mode == "all":
        with tempfile.TemporaryDirectory() as tmpd:
            stats_out = Path(tmpd) / "stats.txt"
            try:
                subprocess.run(
                    ["samtools", "stats", "-@", str(stats_threads), "-F", "0x904", "-o", str(stats_out), str(bam)],
                    check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=stats_timeout
                )
            except subprocess.TimeoutExpired:
                logger.warning("[align_qc] samtools stats timed out; will derive MAPQ/length from sampled reads")
                need_sample_for_stats = True
            except subprocess.CalledProcessError as e:
                logger.warning(f"[align_qc] samtools stats failed ({e}); will derive MAPQ/length from sampled reads")
                need_sample_for_stats = True
            if not need_sample_for_stats:
                for line in stats_out.read_text().splitlines():
                    if line.startswith("MAPQ\t"):
                        _, v, c = line.split("\t")
                        # Cap replication to sample_limit to avoid huge Python lists
                        count = min(int(c), sample_limit)
                        mapq_vals.extend([int(v)] * count)
                    elif line.startswith("RL\t"):
                        _, l, c = line.split("\t")
                        count = min(int(c), sample_limit)
                        read_len_vals.extend([int(l)] * count)
                logger.info(
                    f"[align_qc] samtools stats parsed (capped): mapq_values={len(mapq_vals)}, read_len_values={len(read_len_vals)}"
                )
    elif stats_mode == "skip":
        logger.info("[align_qc] Skipping MAPQ/length stats by configuration")
        need_sample_for_stats = True

    # 3. Sample alignments to compute identity and soft-clip stats
    identity_vals: list[float] = []
    softclip_n = 0
    total_sampled = 0
    # Optional: pre-subsample BAM via samtools view to speed up sampling on huge files
    bam_for_sampling = bam
    samp_spec = os.getenv("FTS_QC_BAM_SAMPLE_FRAC", "").strip()
    if samp_spec:
        try:
            # Accept either raw samtools spec (e.g., '42.001') or a pure fraction (e.g., '0.001')
            if "." in samp_spec and not samp_spec.startswith("0."):
                spec = samp_spec
            else:
                frac = samp_spec
                # derive a spec with default seed 42
                spec = f"42{frac[1:]}" if frac.startswith("0.") else f"42.{frac}"
            with tempfile.TemporaryDirectory() as smpdir:
                samp_bam = Path(smpdir) / "sampled.bam"
                cmd = [
                    "samtools", "view", "-@", str(stats_threads), "-F", "0x904",
                    "-s", spec, "-b", str(bam),
                ]
                with open(samp_bam, "wb") as outfh:
                    subprocess.run(cmd, check=True, stdout=outfh, stderr=subprocess.DEVNULL)
                logger.info(f"[align_qc] Using subsampled BAM for QC sampling: spec={spec}, path={samp_bam}")
                bam_for_sampling = samp_bam
                # Perform sampling within this temp context
                with pysam.AlignmentFile(bam_for_sampling, "rb") as bam_f:
                    for aln in iter_primary(bam_f, sample_limit):
                        total_sampled += 1
                        if need_sample_for_stats:
                            try:
                                mapq_vals.append(int(aln.mapping_quality))
                            except Exception:
                                pass
                            try:
                                if aln.query_length:
                                    read_len_vals.append(int(aln.query_length))
                            except Exception:
                                pass
                        if aln.has_tag("NM") and aln.query_length:
                            identity_vals.append(1 - aln.get_tag("NM") / aln.query_length)
                        if any(op == 4 for op, _ in (aln.cigartuples or [])):
                            softclip_n += 1
        except Exception as e:
            logger.warning(f"[align_qc] Subsample BAM failed ({e}); falling back to direct sampling from full BAM")
            bam_for_sampling = bam

    if bam_for_sampling == bam:
        with pysam.AlignmentFile(bam_for_sampling, "rb") as bam_f:
            for aln in iter_primary(bam_f, sample_limit):
                total_sampled += 1
                if need_sample_for_stats:
                    try:
                        mapq_vals.append(int(aln.mapping_quality))
                    except Exception:
                        pass
                    try:
                        if aln.query_length:
                            read_len_vals.append(int(aln.query_length))
                    except Exception:
                        pass
                if aln.has_tag("NM") and aln.query_length:
                    identity_vals.append(1 - aln.get_tag("NM") / aln.query_length)
                if any(op == 4 for op, _ in (aln.cigartuples or [])):
                    softclip_n += 1
    softclip_pct = percent(softclip_n, total_sampled)
    logger.info(
        f"[align_qc] Sampled primaries: {total_sampled}; identity_n={len(identity_vals)}; softclip_n={softclip_n} ({softclip_pct}%)"
    )
    if need_sample_for_stats:
        logger.info(
            f"[align_qc] Sample-derived MAPQ/length: mapq_values={len(mapq_vals)}, read_len_values={len(read_len_vals)}"
        )

    # 4. Count unique splice junctions (can be expensive on large BEDs)
    unique_juncs = None
    if junc_mode == "skip":
        logger.info("[align_qc] Skipping unique junction count by configuration")
    else:
        eff_limit = junc_limit if (junc_mode == "sample" or junc_limit > 0) else None
        if junc_mode == "sample" and not eff_limit:
            eff_limit = 10_000  # default sample size if not provided
        unique_juncs = count_unique_junctions(bed, sample_limit=eff_limit)
        logger.info(f"[align_qc] Unique junctions: {unique_juncs} (limit={eff_limit or 'none'})")

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
            max_workers=motif_workers,
            max_bed_lines=motif_limit,
        )
        total_motifs = sum(motif_counts.values()) if motif_counts else 0
        logger.info(
            f"[align_qc] Motif counting: keys={len(motif_counts)}, total={total_motifs} (max_bed_lines={motif_limit}, workers={motif_workers})"
        )
    except Exception as e:
        motif_counts = {}
        logging.getLogger(__name__).warning(
            "Splice junction motif counting failed: %s", e
        )
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
