#!/usr/bin/env python3
"""
flair_automate/parse_region_metrics_from_gtf.py

Traverse each region‐directory under a base <regions_dir>,
find the single GTF slice there,
build a gffutils DB, extract per‐gene and per‐transcript metrics,
and write out:

    <region_dir>/gene_summary.csv
    <region_dir>/transcript_summary.csv

If both of those CSVs already exist, the region is skipped.
"""

import os
import glob
import logging
import re
from datetime import datetime

import gffutils
import pandas as pd

def build_db(gtf_path):
    """
    Create or open a persistent gffutils SQLite DB for this GTF.
    """
    db_path = gtf_path + ".db"
    if not os.path.exists(db_path):
        logging.info(f"Building gffutils DB at {db_path}…")
        gffutils.create_db(
            gtf_path, dbfn=db_path, force=True, keep_order=True,
            merge_strategy="merge", sort_attribute_values=True,
            disable_infer_genes=True, disable_infer_transcripts=True
        )
    return gffutils.FeatureDB(db_path, keep_order=True)

def parse_region(db, chrom, start, end):
    """
    Given a gffutils DB & bounds, return two DataFrames:
      - gene_df: one row per gene fully in the interval
      - tx_df:   one row per transcript fully in the interval
    """
    # 1) fetch transcripts fully within [start,end]
    txs = list(db.region(
        seqid=chrom, start=start, end=end,
        featuretype="transcript", completely_within=True
    ))

    # 2) drop any transcript whose child features leak outside
    kept = []
    for tx in txs:
        children = list(db.children(tx, featuretype=None, order_by="start"))
        if all(c.start >= start and c.end <= end for c in children):
            kept.append(tx)
        else:
            logging.debug(f"  dropping {tx.id}: features leak outside {chrom}:{start}-{end}")
    txs = kept

    if not txs:
        raise RuntimeError(f"No fully‐contained transcripts in {chrom}:{start}-{end}")

    # 3) build transcript-level DataFrame
    tx_records = []
    for tx in txs:
        exons    = list(db.children(tx, featuretype="exon", order_by="start"))
        lengths  = [e.end - e.start + 1 for e in exons]
        intervals= [(e.start, e.end) for e in exons]
        tx_records.append({
            "transcript_id":      tx.id,
            "gene_id":            tx.attributes.get("gene_id", [""])[0],
            "chrom":              tx.chrom,
            "tx_start":           tx.start,
            "tx_end":             tx.end,
            "tx_length":          tx.end - tx.start + 1,
            "transcript_biotype": tx.attributes.get("transcript_type", [""])[0],
            "exon_count":         len(exons),
            "exon_lengths":       lengths,
            "exon_intervals":     intervals,
        })
    tx_df = pd.DataFrame(tx_records)

    # 4) build gene‐level DataFrame
    gene_records = []
    for gid in tx_df["gene_id"].unique():
        gene_feat = db[gid]
        # skip genes leaking outside
        if gene_feat.start < start or gene_feat.end > end:
            logging.debug(f"  dropping gene {gid}: span outside region")
            continue

        its   = tx_df[tx_df["gene_id"] == gid]
        lens  = its["tx_length"].tolist()
        exon_counts = its["exon_count"].tolist()

        gene_records.append({
            "gene_id":                  gid,
            "gene_name":                gene_feat.attributes.get("gene_name", [""])[0],
            "chrom":                    gene_feat.chrom,
            "gene_start":               gene_feat.start,
            "gene_end":                 gene_feat.end,
            "gene_length":              gene_feat.end - gene_feat.start + 1,
            "gene_biotype":             gene_feat.attributes.get("gene_type", [""])[0],
            "num_isoforms":             len(its),
            "isoform_ids":              its["transcript_id"].tolist(),
            "transcript_lengths":       lens,
            "mean_tx_length":           pd.Series(lens).mean()   if lens else 0,
            "median_tx_length":         pd.Series(lens).median() if lens else 0,
            "num_exons_per_isoform":    exon_counts,
            "mean_exons_per_isoform":   pd.Series(exon_counts).mean()   if exon_counts else 0,
            "median_exons_per_isoform": pd.Series(exon_counts).median() if exon_counts else 0,
        })
    gene_df = pd.DataFrame(gene_records)

    return gene_df, tx_df

def parse_all_regions(regions_dir, dry_run=False):
    """
    Walk every <regions_dir>/<chrom>_<start>_<end> subdirectory,
    skip if both CSVs already exist,
    find exactly one .gtf there,
    parse & write two summary CSVs.
    """
    logging.info("Parsing GTF metrics under %s", regions_dir)

    # match dir names like chr20_3218000_3250000
    pattern = os.path.join(regions_dir, "*_*_*")
    for region_dir in sorted(glob.glob(pattern)):
        tag      = os.path.basename(region_dir)
        gene_out = os.path.join(region_dir, "gene_summary.csv")
        tx_out   = os.path.join(region_dir, "transcript_summary.csv")

        # skip if already done
        if os.path.exists(gene_out) and os.path.exists(tx_out):
            logging.info("  [SKIP] %s (metrics exist)", tag)
            continue

        # parse chrom/start/end from dir name
        parts = tag.split("_")
        if len(parts) != 3:
            logging.warning("  skipping unrecognized dir %s", tag)
            continue
        chrom, start_s, end_s = parts
        try:
            start, end = int(start_s), int(end_s)
        except ValueError:
            logging.warning("  skipping non‐numeric span in %s", tag)
            continue

        # find the one GTF slice
        gtfs = glob.glob(os.path.join(region_dir, "*.gtf"))
        if len(gtfs) != 1:
            logging.error("  %s: found %d GTF(s), skipping", tag, len(gtfs))
            continue

        gtf = gtfs[0]
        logging.info("  → %s: parsing %s", tag, os.path.basename(gtf))
        if dry_run:
            continue

        # build DB, parse, write CSVs
        db      = build_db(gtf)
        gene_df, tx_df = parse_region(db, chrom, start, end)

        gene_df.to_csv(gene_out, index=False)
        tx_df.to_csv(tx_out,   index=False)
        logging.info("    [✓] wrote %s & %s", os.path.basename(gene_out), os.path.basename(tx_out))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Parse region‐specific GTF metrics under a regions directory"
    )
    parser.add_argument(
        "--regions_dir",
        default=os.path.abspath(os.path.join(
            os.path.dirname(__file__),
            os.pardir, os.pardir, "outputs", "regions"
        )),
        help="Base dir containing <chrom>_<start>_<end> subdirs"
    )
    parser.add_argument(
        "--dry_run", action="store_true",
        help="Do not write CSVs, just log actions"
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %message)s")
    parse_all_regions(args.regions_dir, dry_run=args.dry_run)



