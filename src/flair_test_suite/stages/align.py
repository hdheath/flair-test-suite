# src/flair_test_suite/stages/align.py
from __future__ import annotations
import subprocess
from pathlib import Path

from .base import StageBase
from ..core import PathBuilder
from ..core.input_hash import hash_many


class AlignStage(StageBase):
    """Run `flair align` and publish input‑read count for QC."""
    name = "align"

    @staticmethod
    def _count_reads(fp: Path) -> int:
        """Return number of sequences in FASTA/FASTQ file."""
        if fp.suffix.lower() in {".fa", ".fasta"}:
            return sum(1 for ln in fp.open() if ln.startswith(">"))
        else:
            return sum(1 for _ in fp.open()) // 4

    def build_cmd(self) -> list[str]:
        cfg = self.cfg

        # — resolve core inputs —
        root       = getattr(cfg, "input_root", None) or cfg.run.input_root
        data_dir   = getattr(cfg, "data_dir",   None) or cfg.run.data_dir
        reads_file = getattr(cfg, "reads_file", None) or cfg.run.reads_file
        genome_fa  = getattr(cfg, "genome_fa",  None) or cfg.run.genome_fa

        reads  = (Path(root) / data_dir / reads_file).resolve()
        genome = (Path(root) / data_dir / genome_fa).resolve()

        # — count input reads for QC —
        self._n_input_reads = self._count_reads(reads)

        # — collect signature inputs —
        self._hash_inputs = [reads, genome]

        env        = cfg.run.conda_env
        out_prefix = f"{self.sample}_flair"

        # — parse flags —
        raw_flags = cfg.run.stages[0].flags        # list | str | SimpleNamespace
        flag_parts: list[str] = []

        def _switch(k: str):
            flag_parts.append(k if k.startswith("--") else f"--{k}")

        if isinstance(raw_flags, list):
            flag_parts.extend(raw_flags)

        elif isinstance(raw_flags, str):
            flag_parts.extend(raw_flags.split())

        else:  # table / SimpleNamespace
            for key, val in vars(raw_flags).items():
                if val is False:
                    continue
                if val in ("", None, True):
                    _switch(key)
                elif isinstance(val, (int, float)):
                    _switch(key)
                    flag_parts.append(str(val))
                else:
                    p = Path(val)
                    if not p.is_absolute():
                        p = (Path(root) / data_dir / p).resolve()
                    _switch(key)
                    flag_parts.append(str(p))
                    if p.exists():
                        self._hash_inputs.append(p)

        # — store flags for signature —
        self._flags_components = flag_parts

        # — capture tool version once —
        if not hasattr(self, "_tool_version"):
            try:
                raw = subprocess.check_output(
                    ["conda", "run", "-n", env, "flair", "--version"],
                    text=True
                ).strip()
                self._tool_version = raw.splitlines()[-1] if raw else "flair-unknown"
            except subprocess.CalledProcessError:
                self._tool_version = "flair-unknown"

        # — debug —
        print("[DEBUG] align _hash_inputs:", self._hash_inputs)
        print("[DEBUG] align flags       :", " ".join(flag_parts))

        # — final command —
        return [
            "conda", "run", "-n", env,
            "flair", "align",
            "-g", str(genome),
            "-r", str(reads),
            "-o", out_prefix,
            *flag_parts,
        ]

    @property
    def tool_version(self):
        return getattr(self, "_tool_version", "flair-unknown")

    @property
    def flags_str(self) -> str:
        return " ".join(getattr(self, "_flags_components", []))

    @property
    def input_hashes(self) -> list[str]:
        return hash_many(getattr(self, "_hash_inputs", []))

    def expected_outputs(self):
        base = f"{self.sample}_flair"
        return {
            "bam": Path(f"{base}.bam"),
            "bed": Path(f"{base}.bed"),
        }

    # QC is run by the base class via QC_REGISTRY
    def collect_qc(self, pb):
        return {}


