# src/flair_test_suite/qc/qc_utils.py
# -------------------------------------------------
# Shared helpers used by QC collector modules.

from __future__ import annotations
from pathlib import Path
import subprocess
from typing import Iterator

import pysam
from collections import Counter, defaultdict
from multiprocessing import Pool, cpu_count
import itertools

__all__ = [
    "count_lines",
    "percent",
    "iter_primary",
    "count_unique_junctions",
    "SAMPLE_LIMIT",
]

# ------------------------------------------------------------------
# Simple helpers ----------------------------------------------------

def count_lines(path: Path) -> int:
    """Cheap line‑count using `wc -l`. Raises if the command fails."""
    return int(
        subprocess.check_output(["wc", "-l", str(path)], text=True).split()[0]
    )


def percent(numerator: int, denominator: int, digits: int = 2) -> float:
    """Compute 100·numerator/denominator rounded to *digits* decimals."""
    if denominator == 0:
        return 0.0
    return round(100 * numerator / denominator, digits)

# ------------------------------------------------------------------
# Streaming helper for large BAMs ----------------------------------

# Hard‑coded cap on how many primary alignments we sample when
# computing identities & soft‑clip statistics. Feel free to tweak.
SAMPLE_LIMIT: int = 2_000_000
PRIMARY_ONLY_FLAGMASK = 0x904  # secondary | supplementary | unmapped flags


def iter_primary(bam: pysam.AlignmentFile, limit: int | None = SAMPLE_LIMIT) -> Iterator[pysam.AlignedSegment]:
    """Yield at most *limit* primary alignments from *bam*.

    Primary = not (secondary | supplementary | unmapped).
    """
    n = 0
    for aln in bam.fetch(until_eof=True):
        if aln.flag & PRIMARY_ONLY_FLAGMASK:
            continue
        yield aln
        n += 1
        if limit is not None and n >= limit:
            break

# ------------------------------------------------------------------
# Junction counting from BED12 -------------------------------------

def _parse_blocks(block_sizes: str, block_starts: str) -> tuple[list[int], list[int]]:
    sizes = list(map(int, block_sizes.rstrip(",").split(",")))
    starts = list(map(int, block_starts.rstrip(",").split(",")))
    return sizes, starts


def count_unique_junctions(bed_path: Path, sample_limit: int | None = None) -> int:
    """Return number of unique splice junctions in a BED12 file.

    A junction is defined as (chrom, donor, acceptor, strand). Reads with
    fewer than 2 blocks are skipped.
    *sample_limit* can cap how many lines to parse for very large files;
    default is None (full scan).
    """
    seen: set[tuple[str, int, int, str]] = set()
    with open(bed_path) as fh:
        for idx, line in enumerate(fh):
            if line.startswith("#") or not line.strip():
                continue
            if sample_limit is not None and idx >= sample_limit:
                break
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue  # skip if not BED12
            chrom = cols[0]
            chrom_start = int(cols[1])
            strand = cols[5]
            block_count = int(cols[9])
            if block_count < 2:
                continue
            block_sizes = cols[10]
            block_starts = cols[11]
            sizes, starts = _parse_blocks(block_sizes, block_starts)
            for i in range(block_count - 1):
                donor = chrom_start + starts[i] + sizes[i]
                acceptor = chrom_start + starts[i + 1]
                seen.add((chrom, donor, acceptor, strand))
    return len(seen)

def revcomp(s: str) -> str:
    RC = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return s.translate(RC)[::-1]

def bed_reader(path):
    with open(path) as fh:
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

def chunk(lst, n):
    k, m = divmod(len(lst), n)
    return [
        lst[i*k + min(i, m) : (i+1)*k + min(i+1, m)]
        for i in range(n)
    ]

def _sj_worker(args):
    bed_fields_list, fasta_path, don_ex, don_in, acc_in, acc_ex = args
    fasta = pysam.FastaFile(str(fasta_path))
    ctr = Counter()
    for fields in bed_fields_list:
        for chrom, s, e, strand in introns_from_bed(fields):
            if strand == '+':
                donor    = fasta.fetch(chrom, s - don_ex,   s + don_in)
                acceptor = fasta.fetch(chrom, e - acc_in,   e + acc_ex)
            else:
                donor    = revcomp(fasta.fetch(chrom, e - don_in,   e + don_ex))
                acceptor = revcomp(fasta.fetch(chrom, s - acc_ex,   s + acc_ex))
            ref4 = (donor[-2:] + acceptor[:2]).upper()
            ctr[(strand, ref4)] += 1
    fasta.close()
    return ctr

def count_splice_junction_motifs(
    bed_path,
    fasta_path,
    max_workers=None,
    max_bed_lines=100_000,
    don_ex=2, don_in=2, acc_in=2, acc_ex=2
):
    """
    Count 4-mer splice junction motifs in a BED12 file using the reference genome.
    Returns a Counter: (strand, motif) -> count
    """
    if max_workers is None:
        max_workers = 4  # Default to 4 if not set

    # Read BED lines
    if max_bed_lines is None:
        all_fields = list(bed_reader(bed_path))
    else:
        all_fields = list(itertools.islice(bed_reader(bed_path), max_bed_lines))
    if not all_fields:
        raise ValueError("No usable BED lines!")

    # Group by chromosome
    chrom_groups = defaultdict(list)
    for fields in all_fields:
        chrom_groups[fields[0]].append(fields)

    n_workers = min(max_workers, cpu_count())
    # Build chunks
    if len(chrom_groups) > 1:
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
        single = next(iter(chrom_groups.values()))
        chunks = [c for c in chunk(single, n_workers) if c]

    # Run in parallel
    args = [(chunk, fasta_path, don_ex, don_in, acc_in, acc_ex) for chunk in chunks]
    with Pool(processes=len(chunks)) as pool:
        results = pool.map(_sj_worker, args)

    # Merge results
    total = Counter()
    for ctr in results:
        total.update(ctr)
    return total