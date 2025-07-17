from __future__ import annotations
import subprocess
from pathlib import Path

from .base import StageBase
from ..paths import PathBuilder


class AlignStage(StageBase):
    name = "align"

    # ------------------------------------------------------------------ #
    # Build external command
    # ------------------------------------------------------------------ #
    def build_cmd(self) -> list[str]:
        cfg = self.cfg

        # -------- resolve core inputs ----------------------------------
        root       = getattr(cfg, "input_root", None) or cfg.run.input_root
        data_dir   = getattr(cfg, "data_dir",   None) or cfg.run.data_dir
        reads_file = getattr(cfg, "reads_file", None) or cfg.run.reads_file
        genome_fa  = getattr(cfg, "genome_fa",  None) or cfg.run.genome_fa

        reads  = (Path(root) / data_dir / reads_file).resolve()
        genome = (Path(root) / data_dir / genome_fa).resolve()

        # -------- signature input hashes -------------------------------
        self._input_hashes = [
            PathBuilder.sha256(reads),
            PathBuilder.sha256(genome),
        ]

        env        = cfg.run.conda_env
        out_prefix = f"{self.sample}_flair"

        # -------- flexible flag parsing --------------------------------
        flag_parts: list[str] = []
        raw_flags = cfg.run.stages[0].flags  # list | str | SimpleNamespace

        def _add_switch(k: str):
            flag_parts.append(k if k.startswith("--") else f"--{k}")

        if isinstance(raw_flags, list):
            flag_parts.extend(raw_flags)

        elif isinstance(raw_flags, str):
            flag_parts.extend(raw_flags.split())

        else:  # table / SimpleNamespace
            for key, val in vars(raw_flags).items():
                # switches: true / "" / None
                if val is True or val in ("", None):
                    _add_switch(key)
                # numeric value
                elif isinstance(val, (int, float)):
                    _add_switch(key)
                    flag_parts.append(str(val))
                # string -> maybe a file
                else:
                    p = Path(val)
                    if not p.is_absolute():
                        p = (Path(root) / data_dir / p).resolve()
                    _add_switch(key)
                    flag_parts.append(str(p))
                    if p.exists():
                        self._input_hashes.append(PathBuilder.sha256(p))

        # -------- store flat flag string for signature -----------------
        self._flags_str = " ".join(flag_parts)

        # -------- capture FLAIR version --------------------------------
        if not hasattr(self, "_tool_version"):
            try:
                raw = subprocess.check_output(
                    ["conda", "run", "-n", env, "flair", "--version"],
                    text=True,
                ).strip()
                self._tool_version = raw.splitlines()[-1] if raw else "flair-unknown"
            except subprocess.CalledProcessError:
                self._tool_version = "flair-unknown"

        # -------- final command ----------------------------------------
        return [
            "conda", "run", "-n", env,
            "flair", "align",
            "-g", str(genome),
            "-r", str(reads),
            "-o", out_prefix,
            *flag_parts,
        ]

    # ------------------------------------------------------------------ #
    # Signature helpers
    # ------------------------------------------------------------------ #
    @property
    def tool_version(self):
        return getattr(self, "_tool_version", "flair-unknown")

    @property
    def flags_str(self):
        return self._flags_str

    @property
    def input_hashes(self):
        return self._input_hashes

    # ------------------------------------------------------------------ #
    # Expected outputs
    # ------------------------------------------------------------------ #
    def expected_outputs(self):
        p = f"{self.sample}_flair"
        return {"bam": Path(f"{p}.bam"), "bed": Path(f"{p}.bed")}

    # QC handled globally
    def collect_qc(self, pb):
        return {}

