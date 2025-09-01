# -------------------------------------------------
# QC collector for the region_test stage
# -------------------------------------------------
from __future__ import annotations

import csv
import time
import logging
from collections import defaultdict, Counter
from pathlib import Path
from statistics import mean, median
from typing import Any, Dict, List, Tuple, Iterator

from . import register, write_metrics
from .qc_utils import count_lines, parse_gtf_attributes

__all__ = ["collect"]
logger = logging.getLogger(__name__)


def _stream_gtf(gtf: Path) -> Iterator[List[str]]:
    with open(gtf) as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) >= 9:
                yield cols


def _summarise_gtf(gtf: Path) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    tx_len: Dict[str, int] = {}
    tx_exons: Dict[str, List[int]] = defaultdict(list)
    tx_inter: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    tx_biotype: Dict[str, str] = {}
    tx_gene: Dict[str, str] = {}
    gene_span: Dict[str, Tuple[int, int]] = {}
    gene_biotype: Dict[str, str] = {}
    gene_name: Dict[str, str] = {}
    gene_chrom: Dict[str, str] = {}

    for chrom, _src, ftype, s_s, e_s, _score, _strand, _frame, attrs in _stream_gtf(gtf):
        try:
            s_i, e_i = int(s_s), int(e_s)
        except ValueError:
            continue
        at = parse_gtf_attributes(attrs)
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
        tx_rows.append(
            {
                "transcript_id": tid,
                "gene_id": tx_gene.get(tid, ""),
                "tx_length": length,
                "transcript_biotype": tx_biotype.get(tid, ""),
                "exon_count": len(tx_exons[tid]),
                "exon_lengths": repr(tx_exons[tid]),
                "exon_intervals": repr(tx_inter[tid]),
            }
        )

    gene_rows: List[Dict[str, Any]] = []
    for gid, (g_s, g_e) in gene_span.items():
        iso = [t for t, g in tx_gene.items() if g == gid]
        tx_lengths = [tx_len[t] for t in iso]
        ex_counts = [len(tx_exons[t]) for t in iso]
        gene_rows.append(
            {
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
            }
        )

    return gene_rows, tx_rows


@register("region_test")
def collect(
    manifest: Path,
    out_dir: Path,
    runtime_sec: float | None = None,
    **_ignore,
) -> dict:
    import pandas as pd

    t0 = time.time()
    base_dir = out_dir / "qc"
    details_path = base_dir / "region_details.tsv"
    if not details_path.exists():
        raise FileNotFoundError(f"Missing region details: {details_path}")

    regions = pd.read_csv(details_path, sep="\t")
    region_metrics: List[Dict[str, Any]] = []

    for _, reg in regions.iterrows():
        chrom, start, end = reg["chrom"], int(reg["start"]), int(reg["end"])
        region_tag = f"{chrom}_{start}_{end}"

        bed = out_dir / f"{region_tag}.bed"
        gtf = out_dir / f"{region_tag}.gtf"

        gene_count = transcript_count = exon_count = 0
        biotypes = Counter()

        if gtf.exists() and gtf.stat().st_size > 0:
            gene_rows, tx_rows = _summarise_gtf(gtf)
            gene_count = len(gene_rows)
            transcript_count = len(tx_rows)
            for tx in tx_rows:
                bt = tx.get("transcript_biotype", "")
                if bt:
                    biotypes[bt] += 1
            exon_count = sum(tx.get("exon_count", 0) for tx in tx_rows)

            if gene_rows:
                gene_csv = base_dir / f"{region_tag}_gene_summary.csv"
                with open(gene_csv, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=gene_rows[0].keys())
                    w.writeheader()
                    w.writerows(gene_rows)
            if tx_rows:
                tx_csv = base_dir / f"{region_tag}_transcript_summary.csv"
                with open(tx_csv, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=tx_rows[0].keys())
                    w.writeheader()
                    w.writerows(tx_rows)
        else:
            logger.warning("No GTF for region %s", region_tag)

        bed_lines = 0
        if bed.exists():
            with open(bed) as fh:
                for L in fh:
                    if not L.strip() or L.startswith("#"):
                        continue
                    parts = L.split("\t")
                    if parts[0] != chrom:
                        continue
                    try:
                        b = int(parts[1]) + 1
                        e = int(parts[2])
                    except ValueError:
                        continue
                    if b >= start and e <= end:
                        bed_lines += 1

        span_bp = int(reg["span_bp"]) if "span_bp" in reg else (end - start + 1)
        row: Dict[str, Any] = {
            "region_tag": region_tag,
            "chrom": chrom,
            "start": start,
            "end": end,
            "span_bp": span_bp,
            "gene_count": gene_count,
            "transcript_count": transcript_count,
            "exon_count": exon_count,
            "bed_lines": bed_lines,
        }
        for bt, cnt in biotypes.items():
            row[f"biotype_{bt}"] = cnt
        region_metrics.append(row)

    all_biotypes = set(k for row in region_metrics for k in row.keys() if k.startswith("biotype_"))
    base_keys = [
        "region_tag", "chrom", "start", "end", "span_bp",
        "gene_count", "transcript_count", "exon_count", "bed_lines"
    ]
    keys = base_keys + sorted(all_biotypes)

    for row in region_metrics:
        for k in all_biotypes:
            row.setdefault(k, 0)

    metrics_path = base_dir / "region_metrics.tsv"
    with open(metrics_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(region_metrics)

    metrics = {
        "region_count": len(region_metrics),
        "qc_runtime_sec": round(time.time() - t0, 2),
    }
    if runtime_sec is not None:
        metrics["region_test_runtime_sec"] = round(runtime_sec, 2)

    write_metrics(out_dir, "region_test", metrics)
    return metrics

