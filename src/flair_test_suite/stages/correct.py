# src/flair_test_suite/stages/correct.py
# --------------------------------------
# This stage runs `flair correct` on the BED file produced by the align stage.
# It also uses metadata from the align stage (e.g., input-read count) to enable QC metrics.
# Inherits common orchestration, skip/QC/run logic from StageBase.

from __future__ import annotations

import logging
from pathlib import Path  # for filesystem path operations

from .base import StageBase           # base class providing run() and QC logic
from .stage_utils import (
    resolve_path,
    parse_cli_flags,
    get_stage_config
)
from ..lib.input_hash import hash_many

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
        cfg = self.cfg

        # --- retrieve QC metadata from align stage ---
        meta = self.upstreams["align"].metadata
        self._n_input_reads = meta.get("n_input_reads", meta.get("n_total_reads"))
        if "n_input_reads" not in meta:
            logging.warning(
                "Using legacy metadata key 'n_total_reads' for input read count; consider rerunning align to refresh metadata."
            )
        if self._n_input_reads is None:
            raise RuntimeError(
                "align stage metadata lacks n_input_reads; delete the old signature and rerun align."
            )

        # --- locate files produced by align ---
        align_pb = self.upstreams["align"]
        align_bed = align_pb.stage_dir / f"{self.run_id}_flair.bed"
        if not align_bed.exists():
            logging.warning(
                f"Expected align output BED not found: {align_bed}"
            )
        align_sig = align_pb.signature

        # --- resolve genome reference path ---
        raw_reads = getattr(cfg, "reads_file", None) or cfg.run.reads_file
        resolved = self.resolve_stage_inputs({
            "genome": cfg.run.genome_fa,
            "reads": raw_reads,
        })
        genome = resolved["genome"]
        reads = resolved["reads"]
        self._genome_fa_abs = str(genome)  # <-- Ensure this is set for QC

        # --- inputs for signature ---
        self._hash_inputs = [align_bed, genome]

        # --- parse flags for this stage ---
        flag_parts, extra_inputs = self.resolve_stage_flags()
        self._hash_inputs.extend(extra_inputs)

        # --- include upstream signature in signature inputs ---
        self._hash_inputs.append(align_sig)
        self._flags_components = flag_parts

        # --- warn if no extra flags were supplied ---
        if not flag_parts:
            logging.warning(
                "No extra flags configured for correct stage; using defaults."
            )

        # --- construct final command ---
        flair_version = str(cfg.run.version)
        major_version = int(flair_version.split(".")[0])

        cmd = [
            "flair", "correct",
            "-q", str(align_bed),
            "-o", self.run_id,
            *flag_parts,
        ]
        if major_version < 3:
            cmd.extend(["-g", str(genome)])

        return cmd
        
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



