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

        # --- retrieve QC metadata from upstream ---
        upstream = None
        bed_file = None
        upstream_sig = None

        # Prefer slice if present, otherwise align
        if "slice" in self.upstreams:
            upstream = self.upstreams["slice"]
            bed_file = upstream.stage_dir / "combined_region.bed"
            upstream_sig = upstream.signature
            if not bed_file.exists():
                raise RuntimeError(f"Expected slice output BED not found: {bed_file}")
        elif "align" in self.upstreams:
            upstream = self.upstreams["align"]
            bed_file = upstream.stage_dir / f"{self.run_id}_flair.bed"
            upstream_sig = upstream.signature
            # --- retrieve QC metadata if available ---
            meta = upstream.metadata
            self._n_input_reads = meta.get("n_input_reads", meta.get("n_total_reads"))
            if "n_input_reads" not in meta:
                logging.warning(
                    "Using legacy metadata key 'n_total_reads' for input read count; consider rerunning upstream to refresh metadata."
                )
            if self._n_input_reads is None:
                raise RuntimeError(
                    "Upstream stage metadata lacks n_input_reads; delete the old signature and rerun upstream."
                )

            if not bed_file.exists():
                raise RuntimeError(f"Expected align output BED not found: {bed_file}")
        else:
            raise RuntimeError("correct stage requires upstream 'align' or 'slice' with an existing output BED.")

        # --- resolve genome reference path ---
        raw_reads = getattr(cfg, "reads_file", None) or cfg.run.reads_file
        resolved = self.resolve_stage_inputs({
            "genome": cfg.run.genome_fa,
            "reads": raw_reads,
        })
        genome = resolved["genome"]
        reads = resolved["reads"]
        self._genome_fa_abs = str(genome)

        # --- inputs for signature ---
        self._hash_inputs = [bed_file, genome]

        # --- parse flags for this stage ---
        flag_parts, extra_inputs = self.resolve_stage_flags()
        self._hash_inputs.extend(extra_inputs)

        # --- include upstream signature in signature inputs ---
        self._hash_inputs.append(upstream_sig)
        self._flags_components = flag_parts

        if not flag_parts:
            logging.warning(
                "No extra flags configured for correct stage; using defaults."
            )

        # --- construct final command ---
        flair_version = str(cfg.run.version)
        major_version = int(flair_version.split(".")[0])

        cmd = [
            "flair", "correct",
            "-q", str(bed_file),
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



