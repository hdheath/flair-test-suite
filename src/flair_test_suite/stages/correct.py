# src/flair_test_suite/stages/correct.py
from __future__ import annotations
from pathlib import Path
from .base import StageBase
from ..paths import PathBuilder
from .align import AlignStage

class CorrectStage(StageBase):
    name = "correct"
    requires = ("align",)             # enforce presence
    primary_output_key = "corrected"

    def build_cmd(self) -> list[str]:
        cfg = self.cfg


        align_stage = AlignStage(cfg, self.sample, self.work_dir, self.upstreams)
        align_stage.build_cmd()               # fills _flags_str and _input_hashes
        align_sig  = align_stage.signature
        align_pb   = PathBuilder(cfg.run.work_dir, self.sample, "align", align_sig)
        align_bed  = align_pb.stage_dir / f"{self.sample}_flair.bed"



        root      = getattr(cfg, "input_root", None) or cfg.run.input_root
        data_dir  = getattr(cfg, "data_dir",   None) or cfg.run.data_dir
        genome_fa = getattr(cfg, "genome_fa",  None) or cfg.run.genome_fa
        genome    = (Path(root) / data_dir / genome_fa).resolve()


        # signature inputs
        self._input_hashes = [
            PathBuilder.sha256(align_bed),
            PathBuilder.sha256(genome),
            PathBuilder.sha256_str(align_sig),
        ]


        # ------------------------------------------------------------ flexible flags
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
                    # hash immediately so it’s ready *before* signature is read
                    self._input_hashes.append(PathBuilder.sha256(p))

        self._flags_str = " ".join(flag_parts)            # ➊ already used in signature

        # ---------- signature inputs  (move *below* flag parsing!) -------------
        self._input_hashes.extend([                      # ➋ ensure non‑empty
            PathBuilder.sha256(align_bed),
            PathBuilder.sha256(genome),
            PathBuilder.sha256_str(align_sig),
        ])

        # ---------- final command ----------------------------------------------
        return [
            "conda", "run", "-n", env,
            "flair", "correct",
            "-q", str(align_bed),
            "-o", f"{self.sample}",
            *flag_parts,
        ]


    # ------------------------------------------------------------------ #
    # Files flair‑correct actually writes
    # ------------------------------------------------------------------ #
    def expected_outputs(self):
        base = f"{self.sample}"
        return {
            "corrected":   Path(f"{base}_all_corrected.bed"),
            "inconsistent": Path(f"{base}_all_inconsistent.bed"),
        }

    def collect_qc(self, pb):
        return {}          # handled by qc.correct


