# src/flair_test_suite/qc/align_qc.py
# ---------------------------------
# Enhanced QC collector for the 'align' stage.
#  • Removes deprecated retained_pct key
#  • Adds mapping‑quality, read‑identity, read‑length histograms
#  • Counts unique splice junctions
#  • Records soft‑clipped primary‑alignment counts / percentage
#  • Reports *both* align‑stage runtime and QC runtime
#  • Generates three PNGs in the same QC folder
#
# Compatible with both older (<1.14) and newer samtools versions.
# Requires: pysam, samtools, matplotlib (bundled via flair‑test‑suite env)

from __future__ import annotations
from pathlib import Path
import json
import subprocess
import tempfile
import time
from statistics import mean, median
from ..lib.signature import qc_sidecar_path, load_marker


import matplotlib
matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt  # noqa: E402
import pysam  # noqa: E402

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
    n_input_reads: int,
    genome_fa: str,  # Path to genome FASTA for motif counting
    runtime_sec: float | None = None,  # align‑stage runtime comes from StageBase
) -> dict:
    """QC collector for the *align* stage.

    Outputs:
      • <out_dir>/align_qc.tsv  (metrics)
      • PNG histograms for MAPQ, read identity, read length
    """
    qc_start = time.time()

    # ----------------------------- 1. retained reads -----------------------
    bed = bam.with_suffix(".bed")
    retained = count_lines(bed)
    mapped_pct = percent(retained, n_input_reads)

    # ----------------------------- 2. samtools stats -----------------------
    mapq_vals: list[int] = []
    read_len_vals: list[int] = []

    with tempfile.TemporaryDirectory() as tmpd:
        stats_out = Path(tmpd) / "stats.txt"
        try:
            subprocess.run(
                ["samtools", "stats", "-F", "0x904", "-o", str(stats_out), str(bam)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except subprocess.CalledProcessError:
            with open(stats_out, "w") as fh:
                subprocess.run(
                    ["samtools", "stats", "-F", "0x904", str(bam)],
                    check=True,
                    stdout=fh,
                    stderr=subprocess.DEVNULL,
                )
        for line in stats_out.read_text().splitlines():
            if line.startswith("MAPQ\t"):
                _, v, c = line.split("\t")
                mapq_vals.extend([int(v)] * int(c))
            elif line.startswith("RL\t"):
                _, l, c = line.split("\t")
                read_len_vals.extend([int(l)] * int(c))

    # ----------------------------- 3. identity & softclip ------------------
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

    # ----------------------------- 4. unique junctions ---------------------
    unique_juncs = count_unique_junctions(bed)

    # ----------------------------- 5. plots --------------------------------
    png_root = Path(out_dir)
    png_root.mkdir(parents=True, exist_ok=True)

    def _hist(vals, filename, xlabel):
        if not vals:
            return None
        plt.figure()
        plt.hist(vals, bins=50)
        plt.xlabel(xlabel)
        plt.ylabel("count")
        plt.tight_layout()
        outfile = png_root / filename
        plt.savefig(outfile, dpi=150)
        plt.close()
        return outfile.name

    mapq_png = _hist(mapq_vals, "align_mapq_hist.png", "MAPQ")
    id_png   = _hist([round(v*100,1) for v in identity_vals], "align_identity_hist.png", "Read identity (%)")
    len_png  = _hist(read_len_vals, "align_length_hist.png", "Read length (bp)")

    # Count splice junction motifs (4-mers)
    try:
        motif_counts = count_splice_junction_motifs(
            bed_path=bed,
            fasta_path=Path(genome_fa), 
            max_workers=4  # or get from config/env
        )
    except Exception as e:
        motif_counts = {}
        print(f"Warning: Splice junction motif counting failed: {e}")

    # Save motif counts as a separate JSON file
    motif_counts_str = {f"{k[0]}:{k[1]}": v for k, v in motif_counts.items()}
    with open(Path(out_dir) / "splice_site_motifs.json", "w") as fh:
        json.dump(motif_counts_str, fh, indent=2)

    # Save metrics (without motif counts, or keep as summary)
    metrics = {
        "n_input_reads":   n_input_reads,
        "n_retained_bed":  retained,
        "mapped_pct":      mapped_pct,
        "mean_identity":   round(mean(identity_vals)*100, 2) if identity_vals else 0.0,
        "median_identity": round(median(identity_vals)*100, 2) if identity_vals else 0.0,
        "n_softclip":      softclip_n,
        "softclip_pct":    softclip_pct,
        "unique_junctions": unique_juncs,
        "align_runtime_sec":   round(runtime_sec, 2) if runtime_sec else None,
        "qc_runtime_sec":      round(time.time() - qc_start, 2),
    }

    write_metrics(out_dir, "align", metrics)
    with open(Path(out_dir) / "align_plot_manifest.json", "w") as fh:
        json.dump({"mapq": mapq_png, "identity": id_png, "length": len_png}, fh, indent=2)

    return metrics