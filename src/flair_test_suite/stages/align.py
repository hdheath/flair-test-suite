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

from .stage_utils import estimate_read_count, make_flair_cmd


logger = logging.getLogger(__name__)

class AlignStage(StageBase):
    name = "align"

    def build_cmds(self) -> list[list[str]]:
        cfg = self.cfg

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
        flag_parts, extra_inputs = self.resolve_stage_flags()
        self._hash_inputs.extend(extra_inputs)
        self._flags_components = flag_parts

        # --- capture `flair --version` once per run ---
        if not hasattr(self, "_tool_version"):
            try:
                raw = subprocess.check_output(
                    ["conda", "run", "-n", cfg.run.conda_env, "flair", "--version"],
                    text=True
                ).strip()
                self._tool_version = raw.splitlines()[-1] if raw else "flair-unknown"
            except subprocess.CalledProcessError:
                logger.warning("Could not run `flair --version`; using 'flair-unknown'")
                self._tool_version = "flair-unknown"

        if not flag_parts:
            logger.warning("No extra flags configured for align stage; using defaults")

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
