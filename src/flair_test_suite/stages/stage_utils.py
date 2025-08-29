# src/flair_test_suite/stages/stage_utils.py
# -----------------------------------------
# shared helpers for Stage* classes.

from __future__ import annotations
import gzip
import logging
from pathlib import Path
from typing import Iterable, Iterator, Tuple, List


logger = logging.getLogger(__name__)


def count_reads(fp: Path) -> int:
    """Return number of reads in a FASTA or FASTQ file.

    The original implementation relied solely on the filename suffix to
    distinguish FASTA from FASTQ and could not handle gzipped inputs.  This
    caused ``align`` and ``correct`` QC metrics such as *mapped_pct* to be
    incorrect when the reads were provided as ``fastq`` files (or any gzipped
    variant).  We now detect the format from the first character and transparently
    decompress ``*.gz`` files so either format works.
    """

    opener = gzip.open if fp.suffix == ".gz" else open
    with opener(fp, "rt") as fh:
        first = fh.read(1)
        fh.seek(0)
        if first == ">":
            # FASTA: count header lines
            return sum(1 for ln in fh if ln.startswith(">"))
        if first == "@":
            # FASTQ: total lines divided by four
            return sum(1 for _ in fh) // 4
        # Fallback: guess from filename if the first char was unexpected
        if fp.name.lower().endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz")):
            return sum(1 for ln in fh if ln.startswith(">"))
        return sum(1 for _ in fh) // 4


def estimate_read_count(fp: Path, max_reads: int = 50_000) -> Tuple[int, bool]:
    """Estimate total reads by sampling up to ``max_reads`` records.

    Returns a tuple of (count, is_exact) where ``is_exact`` is True when the
    file ended before the sample limit.
    """

    opener = gzip.open if fp.suffix == ".gz" else open
    total_size = fp.stat().st_size
    with opener(fp, "rt") as fh:
        first = fh.read(1)
        fh.seek(0)
        is_fasta = first == ">" or fp.name.lower().endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz"))
        count = 0
        if is_fasta:
            for line in fh:
                if line.startswith(">"):
                    count += 1
                    if count >= max_reads:
                        break
        else:
            while True:
                header = fh.readline()
                if not header:
                    break
                fh.readline(); fh.readline(); fh.readline()
                count += 1
                if count >= max_reads:
                    break

        # Determine compressed bytes consumed
        consumed = fh.fileobj.tell() if hasattr(fh, "fileobj") else fh.tell()
        eof = fh.readline() == ""  # check if already at EOF

    if count < max_reads and eof:
        return count, True
    # extrapolate
    est = int(count * total_size / consumed) if consumed else 0
    return est, False


def resolve_path(raw: str | Path, *, data_dir: Path) -> Path:
    """
    Expand a path relative to <data_dir> unless it is already absolute.
    """
    p = Path(raw)
    return p if p.is_absolute() else (data_dir / p).resolve()


def make_flair_cmd(
    subcmd: str,
    *,
    genome: Path | str | None = None,
    reads: Iterable[Path | str] | Path | str | None = None,
    bed: Path | str | None = None,
    bam: Path | str | None = None,
    out: Path | str | None = None,
    flags: Iterable[str] | None = None,
) -> List[str]:
    """Return a FLAIR command list with common -g/-r/-q/-b/-o options."""

    cmd: List[str] = ["flair", subcmd]

    if genome:
        cmd.extend(["-g", str(genome)])

    if reads:
        if isinstance(reads, (str, Path)):
            cmd.extend(["-r", str(reads)])
        else:
            cmd.append("-r")
            cmd.extend(str(r) for r in reads)

    if bed:
        cmd.extend(["-q", str(bed)])

    if bam:
        cmd.extend(["-b", str(bam)])

    if out:
        cmd.extend(["-o", str(out)])

    if flags:
        cmd.extend(list(flags))

    return cmd

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

    for k, v in flags_block.items():
        if v is False:  # user explicitly disabled
            logger.warning("Flag '%s' is set to False and will be skipped.", k)
            continue
        if v in (None, "", True):
            _push(k)
        elif isinstance(v, (int, float)):
            _push(k, v)
        elif isinstance(v, str):
            p = resolve_path(v, data_dir=data_dir)
            if p.exists():
                _push(k, p)
                extra_inputs.append(p)
            else:
                logger.warning(
                    "Flag '%s' value '%s' does not resolve to an existing file. Treating as option.",
                    k,
                    v,
                )
                _push(k, v)
        else:
            logger.warning(
                "Flag '%s' has an unrecognized type (%s). Treating as option.", k, type(v)
            )
            _push(k, v)

    return flag_parts, extra_inputs


# ── common stage helpers ---------------------------------------

def collect_upstream_pairs(
    stage_name: str,
    upstreams: dict,
    run_id: str,
    expected_suffix: str,
    regionalized_suffix: str,
    file_check: callable[[Path], bool] = lambda p: p.exists() and p.stat().st_size > 0,
) -> tuple[list[tuple[Path, str]], list[Path], str]:
    """Collect inputs from upstream(s).

    Returns:
      pairs         -> [(path, tag)]
      upstream_sigs -> [signatures]
      mode          -> "regionalized" | "standard"
    """

    pairs: list[tuple[Path, str]] = []
    upstream_sigs: list[Path] = []

    if "regionalize" in upstreams:
        mode = "regionalized"
        reg_pb = upstreams["regionalize"]
        upstream_sigs.append(reg_pb.signature)
        # region_details.tsv is stored under qc/regionalize in newer runs; fall back to root for legacy
        details = reg_pb.stage_dir / "qc" / "regionalize" / "region_details.tsv"
        if not details.exists():
            details = reg_pb.stage_dir / "region_details.tsv"
        if not details.exists():
            raise RuntimeError(f"[{stage_name}] region_details.tsv not found: {details}")

        regions = read_region_details(details)
        for chrom, start, end in regions:
            tag = f"{chrom}_{start}_{end}"
            fpath = None
            for pb in upstreams.values():
                candidate = pb.stage_dir / regionalized_suffix.format(
                    chrom=chrom, start=start, end=end
                )
                if file_check(candidate):
                    fpath = candidate
                    if pb.signature not in upstream_sigs:
                        upstream_sigs.append(pb.signature)
                    break
            if not fpath:
                logging.getLogger(stage_name).warning(
                    "Skipping missing/empty input: %s",
                    regionalized_suffix.format(chrom=chrom, start=start, end=end),
                )
                continue
            pairs.append((fpath, tag))
    else:
        mode = "standard"
        fpath = None
        base_pb = None
        for pb in upstreams.values():
            candidate = pb.stage_dir / f"{run_id}_{expected_suffix}"
            if file_check(candidate):
                base_pb = pb
                fpath = candidate
                break
        if not base_pb or not fpath:
            raise RuntimeError(f"[{stage_name}] requires upstream {expected_suffix}")
        upstream_sigs.append(base_pb.signature)
        pairs.append((fpath, run_id))

    if not pairs:
        raise RuntimeError(f"[{stage_name}] No valid upstream inputs discovered.")

    return pairs, upstream_sigs, mode


def isoform_expected_outputs(run_id: str, first_tag: str | None, regionalized: bool) -> dict[str, Path]:
    if regionalized:
        first = first_tag or run_id
        return {
            "isoforms_bed": Path(f"{first}.isoforms.bed"),
            "isoforms_gtf": Path(f"{first}.isoforms.gtf"),
            "isoforms_bed_pattern": Path("{chrom}_{start}_{end}.isoforms.bed"),
        }
    else:
        return {
            "isoforms_bed": Path(f"{run_id}.isoforms.bed"),
            "isoforms_gtf": Path(f"{run_id}.isoforms.gtf"),
        }


def run_ted_qc(stage_name: str, stage_dir: Path, cfg, upstreams) -> dict:
    # Run TED automatically for collapse/transcriptome; gracefully no-op if
    # required data are unavailable.
    try:
        from ..qc.ted import collect as ted_collect
        import json

        out_dir = stage_dir / "qc" / "ted"
        ted_collect(stage_dir, cfg, upstreams=upstreams, out_dir=out_dir)

        regionalized = any(stage_dir.glob("*_*_*.isoforms.bed"))
        if regionalized:
            map_json = out_dir / "transcriptome_browser" / "region_map.json"
            if not map_json.exists():
                raise RuntimeError("browser plot missing")
            mapping = json.loads(map_json.read_text())
            if not mapping or not all(Path(p).exists() for p in mapping.values()):
                raise RuntimeError("browser plot missing")
        tsv_path = out_dir / "TED.tsv"
        return {"TED": {"tsv": str(tsv_path)}}
    except Exception as e:  # pragma: no cover - logging only
        logging.getLogger(stage_name).warning("TED QC failed: %s", e)
        try:
            (stage_dir / "qc" / "ted" / "TED.tsv").unlink()
        except FileNotFoundError:
            pass
        return {}


def run_sqanti_qc(stage_name: str, stage_dir: Path, cfg, upstreams) -> dict:
    """Run SQANTI QC and return metrics if successful."""
    try:
        from ..qc.sqanti import collect as sqanti_collect

        out_dir = stage_dir / "qc" / "sqanti"
        sqanti_collect(stage_dir, cfg, upstreams=upstreams, out_dir=out_dir)
        tsv = out_dir / "sqanti_results.tsv"
        if not tsv.exists():
            raise RuntimeError("sqanti_results.tsv missing")
        return {"SQANTI": {"tsv": str(tsv)}}
    except Exception as e:  # pragma: no cover - logging only
        logging.getLogger(stage_name).warning("SQANTI QC failed: %s", e)
        try:
            (stage_dir / "qc" / "sqanti" / "sqanti_results.tsv").unlink()
        except FileNotFoundError:
            pass
        return {}


def build_flair_cmds(
    subcmd: str,
    pairs: list[tuple[Path, str]],
    genome: Path,
    reads: list[Path] | None,
    run_id: str,
    flag_parts: list[str],
    regionalized: bool,
    *,
    use_bed: bool = False,
    use_bam: bool = False,
) -> list[list[str]]:
    cmds: list[list[str]] = []
    if regionalized:
        for inp, tag in pairs:
            cmds.append(
                make_flair_cmd(
                    subcmd,
                    bed=inp if use_bed else None,
                    bam=inp if use_bam else None,
                    genome=genome,
                    reads=reads,
                    out=tag,
                    flags=flag_parts,
                )
            )
    else:
        inp, _ = pairs[0]
        cmds.append(
            make_flair_cmd(
                subcmd,
                bed=inp if use_bed else None,
                bam=inp if use_bam else None,
                genome=genome,
                reads=reads,
                out=run_id,
                flags=flag_parts,
            )
        )
    return cmds

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
        logger.warning("Unsupported filetype for filtering: %s", filetype)
        return
    kept = 0
    with open(src) as inf, open(out, "w") as outf:
        for ln, L in enumerate(inf, 1):
            if not L.strip() or L.startswith("#"): continue
            parts = L.rstrip("\n").split("\t")
            if len(parts) <= max(start_i, end_i):
                logger.warning("%s line %s: insufficient cols", src.name, ln)
                continue
            chrom = parts[0]
            try:
                s = int(parts[start_i]) + (1 if add1 else 0)
                e = int(parts[end_i])
            except ValueError:
                logger.warning("%s line %s: bad coords", src.name, ln)
                continue
            if lookup.contains(chrom, s, e):
                outf.write(L)
                kept += 1
    if kept == 0:
        logger.warning("%s empty after filtering", out.name)


def read_region_details(tsv: Path) -> List[tuple[str, str, str]]:
    """Parse region_details.tsv into a list of (chrom, start, end)."""
    regions: List[tuple[str, str, str]] = []
    with tsv.open() as fh:
        next(fh, None)  # skip header
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                logger.warning(
                    "%s line %s: fewer than 3 columns", tsv.name, len(regions) + 2
                )
                continue
            chrom, start, end = parts[0], parts[1], parts[2]
            regions.append((chrom, start, end))
    if not regions:
        raise RuntimeError(f"No regions parsed from {tsv}")
    return regions

def get_stage_config(cfg, name):
    """
    Return the stage config object for the given stage name.
    Raises KeyError if not found.
    """
    for st in cfg.run.stages:
        if st.name == name:
            return st
    raise KeyError(f"Stage '{name}' not found in config.")
