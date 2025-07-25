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
from ..lib.signature import compute_signature
from .stage_utils import (
    count_reads,
    resolve_path,
    parse_cli_flags,
)

class AlignStage(StageBase):
    """
    Run `flair align` and record input-read count for QC.
    Inherits common behavior (run(), signature, skip/QC logic) from StageBase.
    """
    name = "align"

    def build_cmd(self) -> list[str]:
        """
        Construct the command to invoke `flair align` using shared helpers.
        Prepares:
        - self._n_input_reads: for QC collector
        - self._hash_inputs: list of Paths whose contents/signature matter
        - self._flags_components: list of CLI flags for signature
        """
        cfg = self.cfg

        # --- resolve input paths, allowing overrides at top-level or under run ---
        root = Path(getattr(cfg, "input_root", None) or cfg.run.input_root)
        data_dir = Path(getattr(cfg, "data_dir", None) or cfg.run.data_dir)

        reads_file = getattr(cfg, "reads_file", None) or cfg.run.reads_file
        genome_fa = getattr(cfg, "genome_fa", None) or cfg.run.genome_fa

        # Absolute paths using shared resolve_path
        reads = resolve_path(reads_file, root=root, data_dir=data_dir)
        genome = resolve_path(genome_fa, root=root, data_dir=data_dir)

        # --- warn on missing or empty inputs ---
        if not reads.exists():
            warnings.warn(f"Reads file not found: {reads}", UserWarning)
        elif reads.stat().st_size == 0:
            warnings.warn(f"Reads file is empty: {reads}", UserWarning)
        if not genome.exists():
            warnings.warn(f"Genome FASTA not found: {genome}", UserWarning)

        # --- count reads for QC and warn if zero ---
        self._n_input_reads = count_reads(reads)
        if self._n_input_reads == 0:
            warnings.warn(f"No reads counted in {reads}", UserWarning)

        # --- inputs that affect the signature ---
        self._hash_inputs = [reads, genome]

        # --- parse flags from the config for this stage ---
        raw_flags = cfg.run.stages[0].flags
        flag_parts, extra_inputs = parse_cli_flags(raw_flags, root=root, data_dir=data_dir)
        # include any additional files flagged for signature
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
                warnings.warn("Could not run `flair --version`; using 'flair-unknown'", UserWarning)
                self._tool_version = "flair-unknown"

        # --- warn if no extra flags were supplied ---
        if not flag_parts:
            warnings.warn("No extra flags configured for align stage; using defaults", UserWarning)

        # --- construct and return the final command list ---
        out_prefix = f"{self.run_id}_flair"
        return [
            "conda", "run", "-n", cfg.run.conda_env,
            "flair", "align",
            "-g", str(genome),
            "-r", str(reads),
            "-o", out_prefix,
            *flag_parts,
        ]

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
