# src/flair_test_suite/stages/align.py
# -----------------------------------
# This stage runs the `flair align` command on input reads and genome,
# captures the number of input reads for QC, and computes a signature
# based on the inputs, flags, and tool version.

from __future__ import annotations

import warnings          # to emit runtime warnings about config or file issues
import subprocess        # to invoke external commands
from pathlib import Path # for filesystem paths

from .base import StageBase       # base class providing orchestration logic
from ..lib.input_hash import hash_many  # to hash input files


from .stage_utils import (
    count_reads,
    resolve_path,
    parse_cli_flags,
    get_stage_config
)

class AlignStage(StageBase):
    """
    Run `flair align` and record input-read count for QC.
    Inherits common behavior (run(), signature, skip/QC logic) from StageBase.
    """
    name = "align"

    def build_cmd(self) -> list[str]:
        cfg = self.cfg

        # --- resolve input paths ---
        raw_reads = getattr(cfg, "reads_file", None) or cfg.run.reads_file
        resolved = self.resolve_stage_inputs({
            "genome": cfg.run.genome_fa,
            "reads": raw_reads,
        })
        genome = resolved["genome"]
        reads = resolved["reads"]
        self._genome_fa_abs = str(genome)

        # --- count reads for QC and warn if zero ---
        self._n_input_reads = count_reads(reads)
        if self._n_input_reads == 0:
            logging.warning(f"No reads counted in {reads}")

        # --- inputs that affect the signature ---
        self._hash_inputs = [reads, genome]

        # --- parse flags and extra inputs using StageBase helper ---
        flag_parts, extra_inputs = self.resolve_stage_flags()
        self._hash_inputs.extend(extra_inputs)
        self._flags_components = flag_parts

        # --- capture `flair --version` once per run to include in signature ---
        if not hasattr(self, "_tool_version"):
            try:
                raw = subprocess.check_output(
                    ["conda", "run", "-n", cfg.run.conda_env, "flair", "--version"],
                    text=True
                ).strip()
                self._tool_version = raw.splitlines()[-1] if raw else "flair-unknown"
            except subprocess.CalledProcessError:
                logging.warning("Could not run `flair --version`; using 'flair-unknown'")
                self._tool_version = "flair-unknown"

        # --- warn if no extra flags were supplied ---
        if not flag_parts:
            logging.warning("No extra flags configured for align stage; using defaults")

        # --- construct and return the final command list ---
        out_prefix = f"{self.run_id}_flair"
        cmd = [
            "flair", "align",
            "-g", str(genome),
            "-r", str(reads),
            "-o", out_prefix,
            *flag_parts,
        ]
        return cmd

    @property
    def tool_version(self) -> str:
        """Return the cached tool version or a default string."""
        return getattr(self, "_tool_version", "flair-unknown")

    @property
    def flags_str(self) -> str:
        """Stringify flags for signature calculation."""
        return " ".join(getattr(self, "_flags_components", []))

    @property
    def input_hashes(self) -> list[str]:
        """Apply hash_many() to the list of input Paths."""
        return hash_many(getattr(self, "_hash_inputs", []))

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

    def build_cmds(self) -> list[list[str]] | None:
        """Return None as AlignStage uses a single command."""
        return None
