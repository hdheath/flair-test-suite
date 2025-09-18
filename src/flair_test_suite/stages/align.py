# src/flair_test_suite/stages/align.py
# -----------------------------------
# This stage runs the `flair align` command on input reads and genome,
# captures the number of input reads for QC, and computes a signature
# based on the inputs, flags, and tool version.

from __future__ import annotations

import subprocess        # to invoke external commands
import logging
from pathlib import Path # for filesystem paths

from .base import StageBase       # base class providing orchestration logic

from .stage_utils import estimate_read_count, make_flair_cmd, get_stage_config, resolve_path


logger = logging.getLogger(__name__)

class AlignStage(StageBase):
    name = "align"

    def build_cmds(self) -> list[list[str]]:
        cfg = self.cfg

        # Check for user-provided outputs to adopt instead of running FLAIR
        stage_cfg = get_stage_config(cfg, self.name)
        raw_flags = getattr(stage_cfg, "flags", None)
        data_dir = Path(cfg.run.data_dir)
        # Extract --skip value robustly (supports comma-separated BAM,BED and spaces)
        adopt_bam = None
        adopt_bed = None
        import re
        skip_arg: str | None = None
        if isinstance(raw_flags, str):
            # Capture text after --skip (with = or whitespace) up to end or comma followed by another flag
            m = re.search(r"--skip(?:\s+|=)\s*([^,]+(?:,[^,][^,]*)?)(?=\s*,\s*--|$)", raw_flags)
            if m:
                skip_arg = m.group(1).strip()
        elif isinstance(raw_flags, list):
            # Find token that contains --skip, include subsequent non-flag token if present
            for i, tok in enumerate(raw_flags):
                if "--skip" in str(tok):
                    try:
                        parts = re.split(r"(?:\s+|=)", str(tok), maxsplit=1)
                        val = parts[1].strip() if len(parts) > 1 else ""
                    except Exception:
                        val = ""
                    # If next token looks like a path (no leading '-'), append with comma
                    if (i + 1) < len(raw_flags) and not str(raw_flags[i + 1]).strip().startswith("-"):
                        if val:
                            val = f"{val},{str(raw_flags[i + 1]).strip()}"
                        else:
                            val = str(raw_flags[i + 1]).strip()
                    skip_arg = val.strip()
                    break
        # Resolve paths from skip_arg if available
        if skip_arg:
            parts = [p.strip() for p in skip_arg.split(",") if p.strip()]
            if len(parts) == 2:
                adopt_bam = resolve_path(parts[0], data_dir=data_dir)
                adopt_bed = resolve_path(parts[1], data_dir=data_dir)
                logger.info(f"[align] --skip detected; adopting BAM={adopt_bam}, BED={adopt_bed}")

        # --- resolve input paths ---
        raw_reads = getattr(cfg, "reads_file", None) or cfg.run.reads_file
        genome = self.resolve_stage_inputs({"genome": cfg.run.genome_fa})["genome"]

        # normalize reads to list
        if isinstance(raw_reads, (str, Path)):
            raw_reads = [raw_reads]
        if not isinstance(raw_reads, list):
            raise TypeError("reads_file must be a string or list of strings/paths")

        # resolve each read path against data_dir
        resolved_reads = [
            self.resolve_stage_inputs({"reads": r})["reads"] for r in raw_reads
        ]
        self._genome_fa_abs = str(genome)

        # --- estimate reads across all files ---
        total = 0
        all_exact = True
        for r in resolved_reads:
            cnt, exact = estimate_read_count(r)
            total += cnt
            all_exact &= exact
        self._n_input_reads = total
        self._read_count_method = "exact" if all_exact else "estimated"
        if self._n_input_reads == 0:
            logger.warning("No reads counted in any input files: %s", resolved_reads)

        # --- inputs that affect the signature ---
        self._hash_inputs = resolved_reads + [genome]

        # --- parse flags and extra inputs ---
        # Disallow user-provided core IO flags; harness sets these
        reserved = ("r", "reads", "g", "genome", "q", "bed", "b", "bam", "o", "out", "skip")
        flag_parts, extra_inputs = self.resolve_stage_flags(reserved=reserved)
        self._hash_inputs.extend(extra_inputs)
        self._flags_components = flag_parts

        # --- use the version supplied in the config (preferred) ---
        if not hasattr(self, "_tool_version"):
            cfg_ver = getattr(cfg.run, "version", None)
            if cfg_ver:
                # Accept either '2.1.1' or 'flair 2.1.1' from config; normalize sensibly
                s = str(cfg_ver).strip()
                if s.lower().startswith("flair"):
                    self._tool_version = s
                else:
                    self._tool_version = f"flair {s}"
            else:
                logger.warning("No version provided in config.run.version; using 'flair-unknown'")
                self._tool_version = "flair-unknown"

        if not flag_parts:
            logger.warning("No extra flags configured for align stage; using defaults")

        # If both BAM and BED are provided, adopt them and skip running FLAIR.
        # Materialize symlinks in the stage dir so downstreams discover them.
        if adopt_bam and adopt_bed:
            out_prefix = f"{self.run_id}_flair"
            ln_cmd = [
                "bash", "-lc",
                (
                    f"ln -sfn '{adopt_bam}' '{out_prefix}.bam'; "
                    f"ln -sfn '{adopt_bed}' '{out_prefix}.bed'"
                ),
            ]
            logger.info(f"[align] Linking adopted outputs into stage dir with prefix {out_prefix}")
            # Signature inputs reflect overrides + genome + reads + any extra inputs
            self._hash_inputs = [genome, *resolved_reads]
            if adopt_bam:
                self._hash_inputs.append(Path(adopt_bam))
            if adopt_bed:
                self._hash_inputs.append(Path(adopt_bed))
            # No extra flags forwarded; keep the normalized (non-reserved) ones for signature
            reserved = ("r", "reads", "g", "genome", "q", "bed", "b", "bam", "o", "out", "skip")
            flag_parts, extra_inputs = self.resolve_stage_flags(reserved=reserved)
            self._flags_components = flag_parts
            self._hash_inputs.extend(extra_inputs)
            return [ln_cmd]

        # --- construct and return the final command list ---
        out_prefix = f"{self.run_id}_flair"

        # Flair expects one -r arg with comma-separated files
        reads_arg = ",".join(str(r) for r in resolved_reads)

        cmd = make_flair_cmd(
            "align",
            genome=genome,
            reads=reads_arg,
            out=out_prefix,
            flags=flag_parts,
        )
        return [cmd]


    @property
    def tool_version(self) -> str:
        """Return the cached tool version or a default string."""
        return getattr(self, "_tool_version", "flair-unknown")

    def expected_outputs(self) -> dict[str, Path]:
        """
        Map logical output names to file paths within the stage directory.
        'bam' and 'bed' are the primary FLAIR outputs.
        """
        base = f"{self.run_id}_flair"
        return {
            "bam": Path(f"{base}.bam"),
            "bed": Path(f"{base}.bed"),
        }
        
    # QC is invoked automatically by StageBase if a collector is registered
    def collect_qc(self, pb):
        return {}  # no-op here; actual QC logic lives in qc/align_qc.py
