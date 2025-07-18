# src/flair_test_suite/stages/correct.py
from __future__ import annotations
from pathlib import Path

from .base import StageBase
from ..paths import PathBuilder


class CorrectStage(StageBase):
    name = "correct"
    primary_output_key = "corrected"          # <-- QC should use this file

    # ------------------------------------------------------------------ #
    def build_cmd(self) -> list[str]:
        cfg = self.cfg

        # locate align BED
        align_sig = next(
            (Path(cfg.run.work_dir) / self.sample / "align").resolve().iterdir()
        )
        align_bed = (align_sig / f"{self.sample}_flair.bed").resolve()

        # genome FASTA
        root      = (getattr(cfg, "input_root", None) or cfg.run.input_root)
        data_dir  = getattr(cfg, "data_dir", None) or cfg.run.data_dir
        genome_fa = getattr(cfg, "genome_fa", None) or cfg.run.genome_fa
        genome    = (Path(root) / data_dir / genome_fa).resolve()

        self._input_hashes = [
            PathBuilder.sha256(align_bed),
            PathBuilder.sha256(genome),
        ]

        # flexible flags
        env   = cfg.run.conda_env
        flags = cfg.run.stages[1].flags
        flag_parts: list[str] = []

        def _switch(k: str): flag_parts.append(f"--{k}")

        root = Path(root).resolve()
        for k, v in vars(flags).items():
            if v in ("", None, True):
                _switch(k)
            elif v is False:
                continue
            elif isinstance(v, (int, float)):
                _switch(k); flag_parts.append(str(v))
            else:
                p = (root / data_dir / v).resolve()
                _switch(k); flag_parts.append(str(p))
                if p.exists():
                    self._input_hashes.append(PathBuilder.sha256(p))

        self._flags_str = " ".join(flag_parts)

        return [
            "conda", "run", "-n", env,
            "flair", "correct",
            "-q", str(align_bed),
            "-g", str(genome),
            "-o", f"{self.sample}_flair_corr",      # flair appends “_all_*”
            *flag_parts,
        ]

    # ------------------------------------------------------------------ #
    # Files flair‑correct actually writes
    # ------------------------------------------------------------------ #
    def expected_outputs(self):
        base = f"{self.sample}_flair_corr_all"
        return {
            "corrected":   Path(f"{base}_corrected.bed"),
            "inconsistent": Path(f"{base}_inconsistent.bed"),
        }

    def collect_qc(self, pb):
        return {}          # handled by qc.correct


