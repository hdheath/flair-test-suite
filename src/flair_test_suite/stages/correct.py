# src/flair_test_suite/stages/correct.py
# --------------------------------------
# This stage runs `flair correct` on the BED file produced by the align stage.
# It also uses metadata from the align stage (e.g., input-read count) to enable QC metrics.
# Inherits common orchestration, skip/QC/run logic from StageBase.

from __future__ import annotations

import warnings          # to emit runtime warnings
from pathlib import Path  # for filesystem path operations

from .base import StageBase           # base class providing run() and QC logic
from .stage_utils import (
    resolve_path,
    parse_cli_flags,
)

class CorrectStage(StageBase):
    """
    Correct stage for FLAIR pipeline:
    - Takes the BED from align (upstream dependency)
    - Runs `flair correct` with junction and GTF flags
    - Skips or regenerates QC based on existing outputs
    """
    name = "correct"
    requires = ("align",)
    primary_output_key = "corrected"

    def build_cmd(self) -> list[str]:
        """
        Assemble the command to run `flair correct`.
        Also prepares:
        - self._n_input_reads: from align metadata for QC
        - self._hash_inputs: inputs affecting signature (BED + genome + any extra flag files + upstream signature)
        - self._flags_components: CLI flags for signature
        """
        cfg = self.cfg

        # --- retrieve QC metadata from align stage ---
        meta = self.upstreams["align"].metadata
        self._n_input_reads = meta.get("n_input_reads", meta.get("n_total_reads"))
        if "n_input_reads" not in meta:
            warnings.warn(
                "Using legacy metadata key 'n_total_reads' for input read count; consider rerunning align to refresh metadata.",
                UserWarning
            )
        if self._n_input_reads is None:
            raise RuntimeError(
                "align stage metadata lacks n_input_reads; delete the old signature and rerun align."
            )

        # --- locate files produced by align ---
        align_pb = self.upstreams["align"]
        align_bed = align_pb.stage_dir / f"{self.run_id}_flair.bed"
        if not align_bed.exists():
            warnings.warn(
                f"Expected align output BED not found: {align_bed}",
                UserWarning
            )
        align_sig = align_pb.signature

        # --- resolve genome reference path ---
        root = Path(cfg.run.input_root)
        data_dir = Path(cfg.run.data_dir)
        genome = resolve_path(cfg.run.genome_fa, root=root, data_dir=data_dir)
        if not genome.exists():
            warnings.warn(
                f"Genome FASTA not found: {genome}",
                UserWarning
            )

        # --- inputs for signature ---
        self._hash_inputs = [align_bed, genome]

        # --- parse flags for this stage ---
        raw_flags = cfg.run.stages[1].flags
        flag_parts, extra_inputs = parse_cli_flags(raw_flags, root=root, data_dir=data_dir)
        self._hash_inputs.extend(extra_inputs)

        # --- include upstream signature in signature inputs ---
        self._hash_inputs.append(align_sig)
        self._flags_components = flag_parts

        # --- warn if no extra flags were supplied ---
        if not flag_parts:
            warnings.warn(
                "No extra flags configured for correct stage; using defaults.",
                UserWarning
            )

        # --- construct final command ---
        env = cfg.run.conda_env
        return [
            "conda", "run", "-n", env,
            "flair", "correct",
            "-q", str(align_bed),
            "-o", self.run_id,
            *flag_parts,
        ]

    def expected_outputs(self) -> dict[str, Path]:
        """
        Define the output files produced by this stage:
        - <run_id>_all_corrected.bed
        - <run_id>_all_inconsistent.bed
        """
        base = f"{self.run_id}_all"
        return {
            "corrected":    Path(f"{base}_corrected.bed"),
            "inconsistent": Path(f"{base}_inconsistent.bed"),
        }

    def collect_qc(self, pb):
        # no-op here; QC for correct lives in qc/correct_qc.py
        return {}


