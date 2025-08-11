# src/flair_test_suite/qc/slice_qc.py
# -------------------------------------------------
# QC collector for the *slice* stage (combined SliceStage + region-level metrics and plots).
# Requirements:
#  • Python standard library: csv, time, warnings, statistics, collections
#  • pathlib, typing
#  • qc_utils: count_lines
#  • pipeline QC helpers: register, write_metrics
#  • Optional for plots: pandas, numpy, plotly, scipy (install via mamba)
# Summary:
#   Processes the slice manifest and combined_region.gtf/bed outputs to generate:
#     • slice_region_details.tsv (chromosomal regions with span)
#     • gene_summary.csv and transcript_summary.csv (GTF summaries)
#     • region_metrics.tsv (per-region aggregated metrics)
#     • regions_parallel_coordinate.html/png (parallel coords plot)
#   Metrics computed:
#     - Region count, gene count, transcript count
#     - Per-region: gene count, isoform counts, exon stats, biotype proportions, entropy
#     - QC and slice runtimes
# Functions:
#   collect(manifest, out_dir, runtime_sec)
#       Main entry: orchestrates all QC steps for slice stage.
#   _parse_attrs(attr_col)
#       Parse GTF attribute column into dict.
#   _stream_gtf(gtf)
#       Yield non-header GTF fields as lists.
#   _summarise_gtf(gtf)
#       Generate summaries of genes and transcripts from combined_region.gtf.
#   _compute_region_metrics(regions_tsv, gene_csv, tx_csv, max_iso, pc_thresh)
#       Build per-region metric rows and define numeric columns for plotting.
#   _count_lines_for_region(bed, chrom, start, end, col_s, col_e, zero_based)
#       Count features in a BED/JXN file fully inside a region.
#   _try_plot_parcoords(out_dir, region_metrics_tsv, numeric_cols)
#       Produce a Plotly parallel-coordinates chart of region metrics.

from __future__ import annotations
import csv
import time
import warnings
from collections import defaultdict, Counter
from pathlib import Path
from statistics import mean, median
from typing import Any, Dict, List, Tuple

from . import register, write_metrics  # QC registry and output helper
from .qc_utils import count_lines  # simple line counting utility

__all__ = ["collect"]

# ───────────────────────── GTF helpers ────────────────────────────────

def _parse_attrs(attr_col: str) -> Dict[str, str]:
    """
    Convert a GTF attribute column into a key->value dict.
    """
    out: Dict[str, str] = {}
    for chunk in attr_col.rstrip(";").split(";"):
        chunk = chunk.strip()
        if chunk and " " in chunk:
            k, v = chunk.split(" ", 1)
            out[k] = v.strip('"')
    return out


def _stream_gtf(gtf: Path) -> Iterator[List[str]]:
    """
    Yield non-comment lines of a GTF file as lists of fields.
    """
    with open(gtf) as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) >= 9:
                yield cols

# ───────────────────────── summariser ────────────────────────────────

def _summarise_gtf(gtf: Path) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Parse combined_region.gtf to produce gene and transcript summary rows.
    Returns (gene_rows, tx_rows).
    """
    tx_len: Dict[str, int] = {}
    tx_exons: Dict[str, List[int]] = defaultdict(list)
    tx_inter: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
    tx_biotype: Dict[str, str] = {}
    tx_gene: Dict[str, str] = {}
    gene_span: Dict[str, Tuple[int,int]] = {}
    gene_biotype: Dict[str, str] = {}
    gene_name: Dict[str, str] = {}
    gene_chrom: Dict[str, str] = {}

    for chrom, _, ftype, s_s, e_s, *_rest, attrs in _stream_gtf(gtf):
        try:
            s_i, e_i = int(s_s), int(e_s)
        except ValueError:
            continue
        at = _parse_attrs(attrs)
        tid = at.get("transcript_id")
        gid = at.get("gene_id")

        if ftype == "transcript" and tid:
            tx_len[tid] = e_i - s_i + 1
            tx_biotype[tid] = at.get("transcript_type", "")
            if gid:
                tx_gene[tid] = gid
        elif ftype == "exon" and tid:
            tx_exons[tid].append(e_i - s_i + 1)
            tx_inter[tid].append((s_i, e_i))
        elif ftype == "gene" and gid:
            gene_span[gid] = (s_i, e_i)
            gene_biotype[gid] = at.get("gene_type", "")
            gene_name[gid] = at.get("gene_name", "")
            gene_chrom[gid] = chrom

    tx_rows: List[Dict[str, Any]] = []
    for tid, length in tx_len.items():
        tx_rows.append({
            "transcript_id": tid,
            "gene_id": tx_gene.get(tid, ""),
            "tx_length": length,
            "transcript_biotype": tx_biotype.get(tid, ""),
            "exon_count": len(tx_exons[tid]),
            "exon_lengths": repr(tx_exons[tid]),
            "exon_intervals": repr(tx_inter[tid]),
        })

    gene_rows: List[Dict[str, Any]] = []
    for gid, (g_s, g_e) in gene_span.items():
        iso = [t for t, g in tx_gene.items() if g == gid]
        tx_lengths = [tx_len[t] for t in iso]
        ex_counts = [len(tx_exons[t]) for t in iso]
        gene_rows.append({
            "gene_id": gid,
            "gene_name": gene_name.get(gid, ""),
            "chrom": gene_chrom.get(gid, ""),
            "gene_start": g_s,
            "gene_end": g_e,
            "gene_length": g_e - g_s + 1,
            "gene_biotype": gene_biotype.get(gid, ""),
            "num_isoforms": len(iso),
            "isoform_ids": repr(iso),
            "transcript_lengths": repr(tx_lengths),
            "mean_tx_length": mean(tx_lengths) if tx_lengths else 0,
            "median_tx_length": median(tx_lengths) if tx_lengths else 0,
            "num_exons_per_isoform": repr(ex_counts),
            "mean_exons_per_isoform": mean(ex_counts) if ex_counts else 0,
            "median_exons_per_isoform": median(ex_counts) if ex_counts else 0,
        })

    return gene_rows, tx_rows

# ─────────────────── region metrics ─────────────────────────────────

def _compute_region_metrics(
    regions_tsv: Path,
    gene_csv: Path,
    tx_csv: Path,
    max_iso: int = 50,
    pc_thresh: float = 0.0
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """
    Generate per-region metrics and list of numeric columns for plotting.
    """
    import pandas as pd
    from scipy.stats import entropy
    import numpy as np

    regions_df = pd.read_csv(regions_tsv, sep="\t")
    g_df = pd.read_csv(gene_csv, converters={
        "isoform_ids": eval,
        "num_exons_per_isoform": eval,
        "transcript_lengths": eval,
    })
    tx_df = pd.read_csv(tx_csv, converters={
        "exon_lengths": eval,
        "exon_intervals": eval,
    })
    tx_df["transcript_biotype"] = (
        tx_df["transcript_biotype"].astype(str)
        .str.strip().str.lower()
    )

    # Identify top 5 biotypes
    top5 = list(tx_df["transcript_biotype"].value_counts().nlargest(5).index)
    tx2bio = tx_df.set_index("transcript_id")["transcript_biotype"]

    region_rows: List[Dict[str, Any]] = []
    for _, reg in regions_df.iterrows():
        chrom, s, e = reg["chrom"], reg["start"], reg["end"]
        genes = g_df.query("(chrom == @chrom) and (gene_start >= @s) and (gene_end <= @e)")

        if genes.empty:
            # Append zeros for empty regions
            base = {"region_id": f"{chrom}:{s}-{e}", "gene_count": 0, "total_isoforms": 0,
                    "med_isoforms_per_gene": 0, "avg_gene_length": 0,
                    "med_mean_exons_per_isoform": 0, "isoform_entropy": 0,
                    "gene_pc_fraction": 0}
            pct = {f"percent_{bt}": 0.0 for bt in top5}
            base.update(pct); base["percent_other"] = 0.0
            region_rows.append(base)
            continue

        iso_lists = genes["isoform_ids"].tolist()
        iso_ids = [tid for sub in iso_lists for tid in sub]
        if not iso_ids or len(iso_ids) > max_iso:
            continue

        biotypes = [tx2bio.get(t) for t in iso_ids if t in tx2bio.index]
        if not biotypes:
            continue
        cnts = Counter(biotypes); n = sum(cnts.values())
        pct = {f"percent_{bt}": 100 * cnts.get(bt,0)/n for bt in top5}
        pct["percent_other"] = 100 * (n - sum(cnts.get(bt,0) for bt in top5))/n

        pc_frac = genes["gene_biotype"].eq("protein_coding").mean()
        if pc_frac < pc_thresh:
            continue

        iso_counts = genes["num_isoforms"].astype(float).values
        ent = float(entropy(iso_counts/iso_counts.sum(), base=2)) if iso_counts.sum() else 0.0

        region_rows.append({
            "region_id": f"{chrom}:{s}-{e}",
            "gene_count": len(genes),
            "total_isoforms": len(iso_ids),
            "med_isoforms_per_gene": float(genes["num_isoforms"].median()),
            "avg_gene_length": float(genes["gene_length"].mean()),
            "med_mean_exons_per_isoform": float(genes["mean_exons_per_isoform"].median()),
            "isoform_entropy": ent,
            "gene_pc_fraction": pc_frac,
            **pct
        })

    if not region_rows:
        return [], []

    numeric_cols = [
        "gene_count","total_isoforms","med_isoforms_per_gene",
        "avg_gene_length","med_mean_exons_per_isoform",
        "isoform_entropy","gene_pc_fraction"
    ] + sorted([c for c in region_rows[0] if c.startswith("percent_")])

    return region_rows, numeric_cols

# ─────────────────────── per-region line counts ─────────────────────────

def _count_lines_for_region(
    bed: Path, chrom: str, start: int, end: int,
    col_s: int, col_e: int, zero_based: bool
) -> int:
    """
    Count lines in a BED-like file whose features lie fully inside [start,end].
    """
    if not bed.exists():
        return 0
    cnt = 0
    with open(bed) as fh:
        for L in fh:
            if not L.strip() or L.startswith("#"):
                continue
            parts = L.split("\t")
            if parts[0] != chrom:
                continue
            try:
                b = int(parts[col_s]) + (1 if zero_based else 0)
                e = int(parts[col_e])
            except ValueError:
                continue
            if b >= start and e <= end:
                cnt += 1
    return cnt

# ─────────────────────── parallel-coords plot ─────────────────────────

def _try_plot_parcoords(
    out_dir: Path,
    region_metrics_tsv: Path,
    numeric_cols: List[str]
):
    """
    Attempt to generate a parallel-coordinates plot using Plotly.
    Skips if dependencies are missing or data is empty.
    """
    try:
        import pandas as pd
        import numpy as np
        import plotly.graph_objects as go
    except ModuleNotFoundError as e:
        warnings.warn(f"Skipping parcoords plot: missing {e.name}", UserWarning)
        return

    df = pd.read_csv(region_metrics_tsv, sep="\t")
    if df.empty or not numeric_cols:
        warnings.warn("No region metrics → skipping plot", UserWarning)
        return

    dims = []
    for col in numeric_cols:
        mn, mx = df[col].min(), df[col].max()
        ticks = np.linspace(mn, mx, 5)
        ticks_txt = [
            f"{v:.2f}" if col.startswith("percent_") or "fraction" in col else str(int(round(v)))
            for v in ticks
        ]
        label = ("%"+col.split("percent_")[1]) if col.startswith("percent_") else col
        dims.append(dict(range=[mn, mx], tickvals=ticks.tolist(), ticktext=ticks_txt,
                         label=label, values=df[col]))

    trace = go.Parcoords(
        line=dict(color=df["gene_count"], colorscale="Viridis", showscale=True,
                  colorbar=dict(title="Genes")),
        dimensions=dims
    )
    fig = go.Figure(data=[trace])
    fig.update_layout(margin=dict(l=50, r=50, t=50, b=50))
    fig.write_html(out_dir/"regions_parallel_coordinate.html")
    try:
        fig.write_image(out_dir/"regions_parallel_coordinate.png", width=1600, height=800, scale=2)
    except Exception as e:
        warnings.warn(f"Unable to write PNG: {e}", UserWarning)

# ───────────────────────── main collector ───────────────────────────────

@register("regionalize")
def collect(
    manifest: Path,
    out_dir: Path,
    runtime_sec: float | None = None
) -> dict:
    """
    QC collector for per-region RegionalizeStage outputs.
    Processes each region's BED/GTF and reports metrics per region.
    """
    import pandas as pd
    import csv, time, warnings

    t0 = time.time()
    details_path = out_dir / "region_details.tsv"
    if not details_path.exists():
        raise FileNotFoundError(f"Missing region details: {details_path}")

    regions = pd.read_csv(details_path, sep="\t")
    region_metrics = []

    for _, reg in regions.iterrows():
        chrom, start, end = reg["chrom"], reg["start"], reg["end"]
        region_tag = f"{chrom}_{start}_{end}"

        bed = out_dir / f"{region_tag}.bed"
        gtf = out_dir / f"{region_tag}.gtf"
        fa  = out_dir / f"{region_tag}.fa"

        # Compute metrics for this region
        gene_count = transcript_count = exon_count = 0
        biotypes = Counter()
        if gtf.exists():
            gene_rows, tx_rows = _summarise_gtf(gtf)
            gene_count = len(gene_rows)
            transcript_count = len(tx_rows)
            for tx in tx_rows:
                bt = tx.get("transcript_biotype", "")
                if bt: biotypes[bt] += 1
            exon_count = sum(tx.get("exon_count", 0) for tx in tx_rows)

            # --- Write per-region summaries ---
            gene_csv = out_dir / f"{region_tag}_gene_summary.csv"
            tx_csv = out_dir / f"{region_tag}_transcript_summary.csv"
            if gene_rows:
                with open(gene_csv, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=gene_rows[0].keys())
                    w.writeheader()
                    w.writerows(gene_rows)
            if tx_rows:
                with open(tx_csv, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=tx_rows[0].keys())
                    w.writeheader()
                    w.writerows(tx_rows)
        else:
            warnings.warn(f"No GTF for region {region_tag}", UserWarning)

        # Optionally, count lines in BED
        bed_lines = _count_lines_for_region(bed, chrom, start, end, 1, 2, True) if bed.exists() else 0

        region_metrics.append({
            "region_tag": region_tag,
            "chrom": chrom,
            "start": start,
            "end": end,
            "span_bp": reg["span_bp"],
            "gene_count": gene_count,
            "transcript_count": transcript_count,
            "exon_count": exon_count,
            "bed_lines": bed_lines,
            **{f"biotype_{bt}": biotypes[bt] for bt in biotypes}
        })

    # Collect all biotype keys across regions
    all_biotype_keys = set()
    for row in region_metrics:
        all_biotype_keys.update(k for k in row if k.startswith("biotype_"))

    # Build full fieldnames list
    base_keys = [
        "region_tag", "chrom", "start", "end", "span_bp",
        "gene_count", "transcript_count", "exon_count", "bed_lines"
    ]
    keys = base_keys + sorted(all_biotype_keys)

    # Fill missing biotype columns with 0
    for row in region_metrics:
        for k in all_biotype_keys:
            if k not in row:
                row[k] = 0

    # Write per-region metrics
    metrics_path = out_dir / "region_metrics.tsv"
    with open(metrics_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(region_metrics)

    # Final summary metrics
    metrics = {
        "region_count": len(region_metrics),
        "qc_runtime_sec": round(time.time() - t0, 2),
    }
    if runtime_sec is not None:
        metrics["regionalize_runtime_sec"] = round(runtime_sec, 2)
    write_metrics(out_dir, "regionalize", metrics)
    return metrics
