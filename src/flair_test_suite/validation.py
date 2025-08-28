from __future__ import annotations

from collections import Counter
from typing import Iterable

from .config_schema import Config, StageConfig


def _has_manifest_flag(st: StageConfig) -> bool:
    flags = getattr(st, "flags", {}) or {}
    return "manifest" in flags and str(flags["manifest"]).strip() != ""


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
        if n == "regionalize":
            if "align" not in seen:
                raise ValueError("regionalize must appear after align in the TSV list")
        elif n == "correct":
            # correct must follow align (optionally regionalize) and cannot
            # follow any downstream aggregation/terminal stages
            if "align" not in seen:
                raise ValueError("correct must appear after align in the TSV list")
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
            if "align" not in seen:
                raise ValueError(
                    "transcriptome must appear after align (optionally after regionalize) in the TSV list"
                )
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
