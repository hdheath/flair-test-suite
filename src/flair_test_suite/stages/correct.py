# src/flair_test_suite/stages/correct.py
from __future__ import annotations
from pathlib import Path
from .base import StageBase
from ..core import PathBuilder

class CorrectStage(StageBase):
    """
    Run `flair correct` on the BED produced by the upstream align stage.
    The stage is skipped if the signature folder already contains the
    primary output **and** QC side‑car (handled by StageBase).
    """
    name = "correct"
    requires = ("align",)
    primary_output_key = "corrected"

    # ------------------------------------------------------------------ #
    def build_cmd(self) -> list[str]:
        cfg = self.cfg

        # ---------- get align artifacts from upstreams -----------------
        align_pb: PathBuilder = self.upstreams["align"]
        align_bed  = align_pb.stage_dir / f"{self.sample}_flair.bed"
        align_sig  = align_pb.signature           # already stable

        # ---------- resolve shared genome path -------------------------
        root      = cfg.run.input_root
        data_dir  = cfg.run.data_dir
        genome    = (Path(root) / data_dir / cfg.run.genome_fa).resolve()

        # ---------- start collecting input files for hashing -----------
        self._hash_inputs = [align_bed, genome]

        # ---------- parse flag table -----------------------------------
        flags_cfg = cfg.run.stages[1].flags    # second [[run.stages]]
        flag_parts: list[str] = []

        def _add_flag(k: str, v: str | int | None = None):
            flag_parts.append(f"--{k}")
            if v not in (None, ""):
                flag_parts.append(str(v))

        for k, v in vars(flags_cfg).items():
            if v in (None, "", True):
                _add_flag(k)
            elif v is False:
                continue
            elif isinstance(v, (int, float)):
                _add_flag(k, v)
            else:
                fp = (Path(root) / data_dir / v).resolve()
                _add_flag(k, fp)
                self._hash_inputs.append(fp)

        # include the *align* signature itself so different upstream
        # combinations yield a different correct signature folder
        self._hash_inputs.append(align_sig)

        self._flags_components = flag_parts

        # ---------- command --------------------------------------------
        env = cfg.run.conda_env
        return [
            "conda", "run", "-n", env,
            "flair", "correct",
            "-q", str(align_bed),
            "-o", self.sample,               # flair appends “_all_*”
            *flag_parts,
        ]

    # ------------------------------------------------------------------ #
    def expected_outputs(self):
        base = f"{self.sample}_all"
        return {
            "corrected":    Path(f"{base}_corrected.bed"),
            "inconsistent": Path(f"{base}_inconsistent.bed"),
        }

    def collect_qc(self, pb):
        return {}   # handled by qc.correct plug‑in
