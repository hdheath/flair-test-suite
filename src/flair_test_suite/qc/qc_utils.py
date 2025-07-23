# src/flair_test_suite/qc/qc_utils.py
# -------------------------------------------------
# Shared helpers used by QC collector modules.

from __future__ import annotations
from pathlib import Path
import subprocess
from typing import Iterator

import pysam

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
