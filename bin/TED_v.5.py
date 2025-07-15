#!/usr/bin/env python3

import argparse
import sys
import os
import io
import subprocess
from functools import lru_cache, partial
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
import numpy as np

def parse_args():
    p = argparse.ArgumentParser(
        description="Compute precision (bedtools closest) + recall (vectorized overlap) and save summary to TSV."
    )
    p.add_argument('--i', required=True, help="TSV with header: bed prime5 prime3 ref_prime5 ref_prime3")
    p.add_argument('-w', '--window', type=int, default=50, help="±bp window (default: 50)")
    p.add_argument('--summary', required=True, help="Path to write the summary TSV")
    return p.parse_args()

@lru_cache(maxsize=None)
def read_bed6(path):
    return pd.read_csv(
        path, sep='\t', header=None,
        usecols=[0,1,2,5],
        names=['Chrom','Start','End','Strand'],
        dtype={'Chrom':str,'Start':np.int32,'End':np.int32,'Strand':str},
        comment='#'
    )

def prepare_bed6_sorted(path, tmps):
    trimmed = path + ".trimmed.tmp"
    sorted_p = path + ".sorted.tmp"
    with open(trimmed, 'w') as o:
        subprocess.run(["cut", "-f1-6", path], stdout=o, check=True)
    with open(sorted_p, 'w') as o:
        subprocess.run(["bedtools", "sort", "-i", trimmed], stdout=o, check=True)
    os.remove(trimmed)
    tmps.append(sorted_p)
    return sorted_p

def run_closest(a, b):
    cmd = ["bedtools", "closest", "-a", a, "-b", b, "-s", "-d"]
    print("Running command:", ' '.join(cmd), file=sys.stderr)
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    return pd.read_csv(io.StringIO(res.stdout), sep='	', header=None, comment='#')

def extract_distance_and_peak(df, label, max_dist):
    orig_cols = df.shape[1]
    tx_series = df.iloc[:, 3].astype(str)
    dist_series = pd.to_numeric(df.iloc[:, orig_cols - 1], errors='coerce') \
                    .fillna(max_dist + 1).astype(int)
    peak_series = (
        df.iloc[:, 6].astype(str) + '_' +
        df.iloc[:, 7].astype(str) + '_' +
        df.iloc[:, 8].astype(str)
    )
    out = pd.DataFrame({
        'transcript_id': tx_series,
        f'{label}_dist': dist_series,
        f'{label}_peak': peak_series
    })
    return out.groupby('transcript_id', as_index=False).first()

def vectorized_overlap_counts(bed_df, peaks_df, window):
    consumed = np.zeros(len(peaks_df), dtype=bool)
    for (c, s), chunk in bed_df.groupby(['Chrom', 'Strand'], sort=False):
        mask = (peaks_df['Chrom'] == c) & (peaks_df['Strand'] == s)
        idxs = np.flatnonzero(mask)
        if idxs.size == 0:
            continue
        sub = peaks_df.iloc[idxs]
        starts = sub['Start'].to_numpy()
        ends = sub['End'].to_numpy()
        order = np.argsort(starts)
        starts = starts[order]
        ends = ends[order]
        global_i = idxs[order]
        for _, row in chunk.iterrows():
            lo, hi = row['Start'] - window, row['End'] + window
            r = np.searchsorted(starts, hi, side='right')
            if r == 0:
                continue
            hits = np.nonzero(ends[:r] >= lo)[0]
            if hits.size == 0:
                continue
            matched = global_i[hits]
            new_hits = matched[~consumed[matched]]
            consumed[new_hits] = True
    return consumed.sum(), len(peaks_df)

def worker(row, sorted_support, window):
    tmp_local = []
    bed_path = row['bed']
    bed_sorted = prepare_bed6_sorted(bed_path, tmp_local)

    full_bed_df = pd.read_csv(
        bed_path, sep='\t', header=None, comment='#',
        usecols=[0,1,2,3,5],
        names=['Chrom','Start','End','Name','Strand'],
        dtype={'Chrom':str,'Start':np.int32,'End':np.int32,'Name':str,'Strand':str}
    )
    n_tx = len(full_bed_df)

    gene_names = full_bed_df['Name'].str.extract(r'_([^_]+)$')[0]
    n_genes = gene_names.nunique()
    annotated_gene_mask = ~gene_names.str.startswith('chr', na=False)
    annotated_genes = gene_names[annotated_gene_mask].nunique()
    n_annotated_tx = annotated_gene_mask.sum()

    tx_per_gene = n_tx / n_genes if n_genes else 0.0
    tx_per_annot_gene = n_annotated_tx / annotated_genes if annotated_genes else 0.0

    rec = {
        'Primary BED': bed_path,
        'transcripts': n_tx,
        'genes': n_genes,
        'transcripts_per_gene': tx_per_gene,
        'annotated_genes': annotated_genes,
        'transcripts_per_annotated_gene': tx_per_annot_gene
    }

    bed_df = full_bed_df[['Chrom','Start','End','Strand']]

    for col, label in [
        ('prime5', '5'), ('prime3', '3'),
        ('ref_prime5', 'ref5'), ('ref_prime3', 'ref3')
    ]:
        if col in sorted_support:
            dfc = run_closest(bed_sorted, sorted_support[col])
            ddf = extract_distance_and_peak(dfc, label, window)
            m = (ddf[f'{label}_dist'] <= window).sum()
            rec[f'{label}prime_support'] = m
            rec[f'{label}prime_precision'] = m / n_tx if n_tx else 0.0
            if label in ('5', '3'):
                peaks = read_bed6(row[col])
                c, t = vectorized_overlap_counts(bed_df, peaks, window)
                rec[f'{label}prime_peaks_matched'] = c
                rec[f'{label}prime_total_peaks'] = t
                rec[f'{label}prime_recall'] = c / t if t else 0.0
        else:
            if label in ('5', '3'):
                rec.update({
                    f'{label}prime_support': 0,
                    f'{label}prime_precision': 0.0,
                    f'{label}prime_peaks_matched': 0,
                    f'{label}prime_total_peaks': 0,
                    f'{label}prime_recall': 0.0
                })
            else:
                rec.update({
                    f'{label}prime_support': 0,
                    f'{label}prime_precision': 0.0
                })

    for f in tmp_local:
        try:
            os.remove(f)
        except:
            pass

    return rec

def main():
    args = parse_args()

    df = pd.read_csv(args.i, sep=r'\s+', engine='python', comment='#', header=0)
    if 'bed' not in df.columns:
        sys.exit("ERROR: TSV must have a 'bed' column")

    tmp_support = []
    sorted_support = {}
    for col in ('prime5', 'prime3', 'ref_prime5', 'ref_prime3'):
        if col in df.columns and not df[col].isna().all():
            p = df[col].dropna().iloc[0]
            sorted_support[col] = prepare_bed6_sorted(p, tmp_support)

    rows = df.to_dict('records')
    worker_fn = partial(worker, sorted_support=sorted_support, window=args.window)
    with ProcessPoolExecutor(max_workers=12) as exe:
        records = list(exe.map(worker_fn, rows))

    for f in tmp_support:
        try:
            os.remove(f)
        except:
            pass

    out = pd.DataFrame(records)

    def safe_f1(p, r):
        return 2 * p * r / (p + r) if (p + r) > 0 else 0.0

    out['5prime_f1'] = [safe_f1(p, r) for p, r in zip(out['5prime_precision'], out['5prime_recall'])]
    out['3prime_f1'] = [safe_f1(p, r) for p, r in zip(out['3prime_precision'], out['3prime_recall'])]
    out['mean_f1'] = (out['5prime_f1'] + out['3prime_f1']) / 2

    out['sort_score'] = out['mean_f1']
    out.sort_values(by='sort_score', ascending=False, inplace=True)

    cols = [
        'Primary BED', 'transcripts', 'genes', 'transcripts_per_gene',
        'annotated_genes', 'transcripts_per_annotated_gene',
        '5prime_support','5prime_precision','5prime_peaks_matched','5prime_total_peaks','5prime_recall',
        '3prime_support','3prime_precision','3prime_peaks_matched','3prime_total_peaks','3prime_recall',
        'ref5prime_support','ref5prime_precision',
        'ref3prime_support','ref3prime_precision',
        '5prime_f1', '3prime_f1', 'mean_f1'
    ]

    out[cols].to_csv(args.summary, sep='\t', index=False)
    print(f"Summary written to {args.summary}")

if __name__ == '__main__':
    main()


'''
python TED_v.5.py --i /private/groups/brookslab/hdheath/projects/rotation/testing/yeast_precision_recall.tsv --summary Yeast_aim2-preci-recall.tsv 

python TED_v.5.py --i /private/groups/brookslab/hdheath/projects/rotation/testing/human_precision_recall.tsv --summary Human_aim2-preci-recall.tsv 

python TED_v.5.py --i /private/groups/brookslab/hdheath/projects/rotation/test_suite/aim3-yeast_precision_recall.tsv --summary aim3-Yeast-preci-recall.tsv


'''