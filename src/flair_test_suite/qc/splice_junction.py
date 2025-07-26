#!/usr/bin/env python3
"""
bed_reference_motifs_parallel.py

Extract 4‑mer splice‑site motifs from the reference genome based on a BED12,
parallelized across up to MAX_WORKERS processes, and process at most
MAX_BED_LINES reads from the BED.

Usage:
    python bed_reference_motifs_parallel.py
"""

import sys
import itertools
from pathlib import Path
from collections import Counter, defaultdict
from multiprocessing import Pool, cpu_count
import pysam

# ── user paths ──────────────────────────────────────────────────────
BED_PATH   = Path("/private/groups/brookslab/hdheath/projects/flair-test-suite/flair-test-suite/outputs/WTC11_demo/align/3d23cbcd128f/WTC11_demo_flair.bed")
FASTA_PATH = Path("/private/groups/brookslab/hdheath/projects/flair-test-suite/flair-test-suite/tests/data/GRCh38.primary_assembly.genome.fa")

# 7e2d8d3a6be0
# 3d23cbcd128f

# ── processing caps ─────────────────────────────────────────────────────────
MAX_WORKERS   = 8       # max parallel processes
MAX_BED_LINES = 100_000  # at most this many BED lines (None to disable)

# ── motif flank sizes ─────────────────────────────────────────────────────────
DON_EX, DON_IN  = 2, 2   # donor: exonic/intronic
ACC_IN, ACC_EX  = 2, 2   # acceptor: intronic/exonic

# ── helper ────────────────────────────────────────────────────────────────────
RC = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def revcomp(s: str) -> str:
    return s.translate(RC)[::-1]

# ── parse BED12 introns ───────────────────────────────────────────────────────
def bed_reader(path: Path):
    with path.open() as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("track","browser","#")):
                continue
            yield ln.rstrip().split('\t')

def introns_from_bed(fields):
    chrom      = fields[0]
    chromStart = int(fields[1])
    strand     = fields[5]
    sizes      = list(map(int, fields[10].rstrip(',').split(',')))
    starts     = list(map(int, fields[11].rstrip(',').split(',')))
    for i in range(len(sizes) - 1):
        exon1_end   = chromStart + starts[i]     + sizes[i]
        exon2_start = chromStart + starts[i+1]
        yield chrom, exon1_end, exon2_start, strand

# ── chunk helper ─────────────────────────────────────────────────────────────
def chunk(lst, n):
    k, m = divmod(len(lst), n)
    return [
        lst[i*k + min(i, m) : (i+1)*k + min(i+1, m)]
        for i in range(n)
    ]

# ── worker ──────────────────────────────────────────────────────────────────
def worker(bed_fields_list):
    fasta = pysam.FastaFile(str(FASTA_PATH))
    ctr = Counter()
    for fields in bed_fields_list:
        for chrom, s, e, strand in introns_from_bed(fields):
            if strand == '+':
                donor    = fasta.fetch(chrom, s - DON_EX,   s + DON_IN)
                acceptor = fasta.fetch(chrom, e - ACC_IN,   e + ACC_EX)
            else:
                donor    = revcomp(fasta.fetch(chrom, e - DON_IN,   e + DON_EX))
                acceptor = revcomp(fasta.fetch(chrom, s - ACC_EX,   s + ACC_EX))
            ref4 = (donor[-2:] + acceptor[:2]).upper()
            ctr[(strand, ref4)] += 1
    fasta.close()
    return ctr

# ── main ────────────────────────────────────────────────────────────────────
def main():
    # stream up to MAX_BED_LINES
    if MAX_BED_LINES is None:
        all_fields = list(bed_reader(BED_PATH))
    else:
        all_fields = list(itertools.islice(bed_reader(BED_PATH), MAX_BED_LINES))

    if not all_fields:
        print("No usable BED lines!", file=sys.stderr)
        sys.exit(1)

    # group fields by chromosome
    chrom_groups = defaultdict(list)
    for fields in all_fields:
        chrom_groups[fields[0]].append(fields)

    # determine number of workers
    n_workers = min(MAX_WORKERS, cpu_count())

    # build chunks: by chromosome or split single-chrom list
    if len(chrom_groups) > 1:
        # distribute chromosomes roughly across workers
        chroms = list(chrom_groups.keys())
        chrom_chunks = chunk(chroms, min(n_workers, len(chroms)))
        chunks = []
        for cchunk in chrom_chunks:
            group = []
            for chrom in cchunk:
                group.extend(chrom_groups[chrom])
            if group:
                chunks.append(group)
    else:
        # only one chromosome: split its list into n_workers pieces
        single = next(iter(chrom_groups.values()))
        chunks = [c for c in chunk(single, n_workers) if c]

    # run in parallel
    with Pool(processes=len(chunks)) as pool:
        results = pool.map(worker, chunks)

    # merge results
    total = Counter()
    for ctr in results:
        total.update(ctr)

    # output summary
    print("strand\tref4\tcount")
    for (strand, motif), count in total.most_common():
        print(f"{strand}\t{motif}\t{count}")

if __name__ == "__main__":
    main()
