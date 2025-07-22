# src/flair_test_suite/stages/correct.py
# --------------------------------------
# This stage runs `flair correct` on the BED file produced by the align stage.
# It also uses metadata from the align stage (e.g., input-read count) to enable QC metrics.
# Inherits common orchestration, skip/QC/run logic from StageBase.

from __future__ import annotations
import warnings          # to emit runtime warnings
from pathlib import Path  # for filesystem path operations

from .base import StageBase  # base class providing run() and QC logic
from ..lib import PathBuilder  # constructs stage directories


class CorrectStage(StageBase):
    """
    Correct stage for FLAIR pipeline:
    - Takes the BED from align (upstream dependency)
    - Runs `flair correct` with junction and GTF flags
    - Skips or regenerates QC based on existing outputs
    """
    # Unique identifier used by StageBase and QC registry
    name = "correct"
    # This stage depends on the align stage having completed
    requires = ("align",)
    # Which key in expected_outputs() is the primary file to check
    primary_output_key = "corrected"

    # ------------------------------------------------------------------ #
    def build_cmd(self) -> list[str]:
        """
        Assemble the command to run `flair correct`.
        Also prepares:
        - self._n_input_reads: from align metadata for QC
        - self._hash_inputs: inputs affecting signature (BED + genome)
        - self._flags_components: CLI flags for signature
        """
        cfg = self.cfg

        # --- retrieve QC metadata from align stage ---
        meta = self.upstreams["align"].metadata
        # prefer new key, fallback for old markers
        self._n_input_reads = meta.get("n_input_reads", meta.get("n_total_reads"))
        if "n_input_reads" not in meta:
            warnings.warn(
                "Using legacy metadata key 'n_total_reads' for input read count; consider rerunning align to refresh metadata.",
                UserWarning
            )
        if self._n_input_reads is None:
            # bail out if no metadata present
            raise RuntimeError(
                "align stage metadata lacks n_input_reads; "
                "delete the old signature folder and rerun align."
            )

        # --- locate files produced by align in its stage directory ---
        align_pb: PathBuilder = self.upstreams["align"]
        align_bed = align_pb.stage_dir / f"{self.run_id}_flair.bed"
        if not align_bed.exists():
            warnings.warn(
                f"Expected align output BED not found: {align_bed}",
                UserWarning
            )
        # reuse align signature so signature reflects upstream outputs
        align_sig = align_pb.signature

        # --- resolve genome reference path from config ---
        root = cfg.run.input_root
        data_dir = cfg.run.data_dir
        genome = (Path(root) / data_dir / cfg.run.genome_fa).resolve()
        if not genome.exists():
            warnings.warn(
                f"Genome FASTA not found: {genome}",
                UserWarning
            )

        # --- inputs contributing to this stage's signature ---
        self._hash_inputs = [align_bed, genome]

        # --- parse flags for this stage from the TOML config ---
        flags_cfg = cfg.run.stages[1].flags  # second [[run.stages]] block
        flag_parts: list[str] = []

        def _add_flag(k: str, v: str | int | None = None):
            # always prefix keys with '--'
            flag_parts.append(f"--{k}")
            if v not in (None, ""):
                flag_parts.append(str(v))

        # iterate all key/value pairs in the config block
        for k, v in vars(flags_cfg).items():
            if v in (None, "", True):
                # switch flag with no argument
                _add_flag(k)
            elif v is False:
                # user disabled this flag explicitly
                continue
            elif isinstance(v, (int, float)):
                # numeric flag: add key and value
                _add_flag(k, v)
            else:
                # treat value as a relative path under root/data_dir
                fp = (Path(root) / data_dir / v).resolve()
                _add_flag(k, fp)
                if fp.exists():
                    self._hash_inputs.append(fp)
                else:
                    warnings.warn(
                        f"Flag file {fp} does not exist; skipping it in signature.",
                        UserWarning
                    )

        # include the align signature in this stage's signature
        self._hash_inputs.append(align_sig)
        # store flags for signature computation
        self._flags_components = flag_parts

        # --- warn if no extra flags were supplied ---
        if not flag_parts:
            warnings.warn(
                "No extra flags configured for correct stage; using defaults.",
                UserWarning
            )

        # --- construct final command to execute ---
        env = cfg.run.conda_env
        return [
            "conda", "run", "-n", env,
            "flair", "correct",
            "-q", str(align_bed),       # input BED from align
            "-o", self.run_id,          # prefix for output files
            *flag_parts,                 # include parsed flags
        ]

    # ------------------------------------------------------------------ #
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

# Placeholder: additional stages can be added below
# e.g., class CollapseStage(StageBase): ...

