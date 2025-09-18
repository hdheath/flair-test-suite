from __future__ import annotations

from collections import Counter
from typing import Iterable

from .config_schema import Config, StageConfig


def _has_manifest_flag(st: StageConfig) -> bool:
    flags = getattr(st, "flags", None)
    if isinstance(flags, dict):
        v = flags.get("manifest")
        return bool(v and str(v).strip())
    if isinstance(flags, str):
        toks = [t.strip() for t in flags.split(',') if t.strip()]
    elif isinstance(flags, list):
        toks = [str(t).strip() for t in flags if str(t).strip()]
    else:
        toks = []
    for tok in toks:
        t = tok.lstrip('-')
        if t.startswith('manifest='):
            return True
    return False


def _flags_dict(st: StageConfig) -> dict:
    """Best-effort parse of a stage's flags into a dict-like mapping.

    Supports dict, list[str], or comma-separated string. Values are strings when
    provided as --key=value or "--key value"; bare flags map to True.
    """
    flags = getattr(st, "flags", None)
    if isinstance(flags, dict):
        # normalize keys to bare form
        return {str(k).lstrip('-'): v for k, v in flags.items()}
    toks: list[str]
    if isinstance(flags, str):
        toks = [t.strip() for t in flags.split(',') if t.strip()]
    elif isinstance(flags, list):
        toks = [str(t).strip() for t in flags if str(t).strip()]
    else:
        toks = []
    out: dict[str, object] = {}
    import shlex
    for tok in toks:
        parts = shlex.split(tok)
        if not parts:
            continue
        head = parts[0]
        if '=' in head:
            k, v = head.lstrip('-').split('=', 1)
            out[k] = v
        else:
            k = head.lstrip('-')
            if len(parts) > 1 and not parts[1].startswith('-'):
                out[k] = parts[1]
            else:
                out[k] = True
    return out  # type: ignore[return-value]


def validate_stage_order(cfg: Config) -> None:
    """Validate that stages listed in cfg.run.stages are executable in-order.

    Rules (implicit dependencies derived from stage implementations):
      - regionalize: must follow align
      - correct:     must follow align (may optionally follow regionalize too)
      - collapse:    must follow correct
      - transcriptome: must follow align (regionalize is optional and itself follows align)
      - combine:     must follow collapse or transcriptome, unless a manifest flag is provided
      - quantify:    must follow combine OR transcriptome OR collapse

    Also enforces:
      - No duplicate stage names (duplicates would be collapsed later and are ambiguous)
    """

    stages = list(getattr(cfg.run, "stages", []) or [])
    names = [getattr(s, "name", "") for s in stages]

    # 1) Deduplicate check
    counts = Counter(names)
    dups = [n for n, c in counts.items() if c > 1]
    if dups:
        raise ValueError(
            "Duplicate stages detected in config: "
            + ", ".join(sorted(dups))
            + ". Each stage should appear at most once in the TSV order."
        )

    # 2) Order constraints
    seen: set[str] = set()
    for st in stages:
        n = st.name
        if n == "region_test":
            # Allow starting at region_test when both BAM and BED are supplied via
            # stage flags; otherwise require align first.
            flags = _flags_dict(st)
            have_overrides = bool(flags.get("bam") and flags.get("bed"))
            if not have_overrides and "align" not in seen:
                raise ValueError("region_test must appear after align, or provide both --bam and --bed")
        elif n == "correct":
            # correct must follow align (optionally regionalize) and cannot
            # follow any downstream aggregation/terminal stages
            flags = _flags_dict(st)
            have_bed_override = bool(flags.get("bed"))
            if not have_bed_override and "align" not in seen:
                raise ValueError("correct must appear after align, or provide --bed")
            forbidden = {"transcriptome", "collapse", "combine", "quantify"}
            bad = forbidden.intersection(seen)
            if bad:
                raise ValueError(
                    "correct cannot follow these stages: " + ", ".join(sorted(bad))
                )
        elif n == "collapse":
            if "correct" not in seen:
                raise ValueError("collapse must appear after correct in the TSV list")
        elif n == "transcriptome":
            # transcriptome requires align and optionally regionalize, but must
            # NOT follow correct. Enforce align seen and no prior correct.
            flags = _flags_dict(st)
            have_bam_override = bool(flags.get("bam"))
            if not have_bam_override and "align" not in seen:
                raise ValueError("transcriptome must appear after align (or provide --bam); region_test may precede it")
            if "correct" in seen:
                raise ValueError(
                    "transcriptome cannot follow correct; remove 'correct' or place 'transcriptome' before it"
                )
        elif n == "combine":
            if ("collapse" not in seen and "transcriptome" not in seen) and not _has_manifest_flag(st):
                raise ValueError(
                    "combine must appear after collapse or transcriptome, or provide flags.manifest"
                )
        elif n == "quantify":
            if not ("combine" in seen or "transcriptome" in seen or "collapse" in seen):
                raise ValueError(
                    "quantify must appear after combine, transcriptome, or collapse in the TSV list"
                )

        seen.add(n)
