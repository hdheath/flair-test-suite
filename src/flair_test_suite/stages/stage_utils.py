# src/flair_test_suite/stages/stage_utils.py
# -----------------------------------------
# Hand‑rolled shared helpers for Stage* classes.
# Keeping them here (rather than in core/) lets us import without a circular
# dependency, since StageBase already lives in stages/.

from __future__ import annotations
import warnings
from pathlib import Path
from typing import Iterable, Iterator, Tuple, List

# ── I/O helpers -----------------------------------------------------

def count_reads(fp: Path) -> int:
    """Return # reads in a FASTA/FASTQ file (cheap header scan)."""
    if fp.suffix.lower() in {".fa", ".fasta"}:
        return sum(1 for ln in fp.open() if ln.startswith(">"))
    return sum(1 for _ in fp.open()) // 4


def resolve_path(raw: str | Path, *, root: Path, data_dir: Path) -> Path:
    """
    Expand a path relative to <root>/<data_dir> unless it is already absolute.
    """
    p = Path(raw)
    return p if p.is_absolute() else (root / data_dir / p).resolve()

# ── flag parsing helpers -------------------------------------------

def iter_stage_flags(flags_block) -> Iterator[Tuple[str, object]]:
    """
    Yield (key, value) from any TOML [[run.stages]].flags variant:
      • table  ‑> attrs in declaration order
      • list   ‑> each item becomes (flag, True)
      • string ‑> whitespace‑split tokens
    """
    if flags_block is None:
        return
    if isinstance(flags_block, list):          # ["--foo", "--bar", …]
        for tok in flags_block:
            yield tok.lstrip("-"), True
    elif isinstance(flags_block, str):         # "--foo  --bar=baz"
        for tok in flags_block.split():
            if "=" in tok:
                k, v = tok.split("=", 1)
                yield k.lstrip("-"), v
            else:
                yield tok.lstrip("-"), True
    else:                                      # TOML table
        for k, v in vars(flags_block).items():
            yield k, v


def parse_cli_flags(
    flags_block,
    *,
    root: Path,
    data_dir: Path,
) -> Tuple[List[str], List[Path]]:
    """
    Convert a flags block into:
      • flag_parts   (for subprocess command)
      • extra_inputs (files that should be hashed)
    """
    flag_parts: List[str] = []
    extra_inputs: List[Path] = []

    def _push(k: str, v: str | int | None = None):
        if len(k) == 1:
            flag_parts.append(f"-{k}")
        else:
            flag_parts.append(f"--{k}")
        if v not in (None, "", True):
            flag_parts.append(str(v))

    for k, v in iter_stage_flags(flags_block):
        if v is False:           # user explicitly disabled
            continue
        if v in (None, "", True):
            _push(k)
        elif isinstance(v, (int, float)):
            _push(k, v)
        else:                    # assume path‑like
            p = resolve_path(v, root=root, data_dir=data_dir)
            _push(k, p)
            extra_inputs.append(p)

    return flag_parts, extra_inputs

def filter_file_by_regions(src: Path, out: Path, lookup, filetype: str):
    """
    Filter lines from src by region containment and write to out.
    filetype: "bed", "tab", "gtf"
    """
    if filetype == "bed":
        start_i, end_i, add1 = 1, 2, True
    elif filetype == "tab":
        start_i, end_i, add1 = 1, 2, False
    elif filetype == "gtf":
        start_i, end_i, add1 = 3, 4, False
    else:
        warnings.warn(f"Unsupported filetype for filtering: {filetype}", UserWarning)
        return
    kept = 0
    with open(src) as inf, open(out, "w") as outf:
        for ln, L in enumerate(inf, 1):
            if not L.strip() or L.startswith("#"): continue
            parts = L.rstrip("\n").split("\t")
            if len(parts) <= max(start_i, end_i):
                warnings.warn(f"{src.name} line {ln}: insufficient cols", UserWarning)
                continue
            chrom = parts[0]
            try:
                s = int(parts[start_i]) + (1 if add1 else 0)
                e = int(parts[end_i])
            except ValueError:
                warnings.warn(f"{src.name} line {ln}: bad coords", UserWarning)
                continue
            if lookup.contains(chrom, s, e):
                outf.write(L)
                kept += 1
    if kept == 0:
        warnings.warn(f"{out.name} empty after filtering", UserWarning)