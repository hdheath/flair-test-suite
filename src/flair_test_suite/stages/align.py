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
from ..lib import PathBuilder    # builder for output directories
from ..lib.input_hash import hash_many  # to hash input files

class AlignStage(StageBase):
    """
    Run `flair align` and record input-read count for QC.
    Inherits common behavior (run(), signature, skip/QC logic) from StageBase.
    """
    # Unique name used for registry lookup and QC tagging
    name = "align"

    @staticmethod
    def _count_reads(fp: Path) -> int:
        """
        Return the number of reads in a FASTA/FASTQ file.
        - FASTA: count lines starting with '>'.
        - FASTQ: assume 4 lines per record, divide total lines by 4.
        """
        if fp.suffix.lower() in {".fa", ".fasta"}:
            # FASTA mode: header lines start with '>'
            return sum(1 for ln in fp.open() if ln.startswith(">"))
        else:
            # FASTQ mode: 4 lines per read
            return sum(1 for _ in fp.open()) // 4

    def build_cmd(self) -> list[str]:
        """
        Construct the command to invoke `flair align`.
        Also prepares:
        - self._n_input_reads: for QC collector
        - self._hash_inputs: list of Paths whose contents/signature matter
        - self._flags_components: list of CLI flags for signature
        """
        cfg = self.cfg

        # --- resolve input paths, allowing override at top-level or under run ---
        root       = getattr(cfg, "input_root", None) or cfg.run.input_root
        data_dir   = getattr(cfg, "data_dir",   None) or cfg.run.data_dir
        reads_file = getattr(cfg, "reads_file", None) or cfg.run.reads_file
        genome_fa  = getattr(cfg, "genome_fa",  None) or cfg.run.genome_fa

        # Absolute paths to the reads and genome files
        reads  = (Path(root) / data_dir / reads_file).resolve()
        genome = (Path(root) / data_dir / genome_fa).resolve()

        # --- warn on missing or empty inputs ---
        if not reads.exists():
            warnings.warn(f"Reads file not found: {reads}", UserWarning)
        elif reads.stat().st_size == 0:
            warnings.warn(f"Reads file is empty: {reads}", UserWarning)
        if not genome.exists():
            warnings.warn(f"Genome FASTA not found: {genome}", UserWarning)

        # --- count reads for QC and warn if zero ---
        self._n_input_reads = self._count_reads(reads)
        if self._n_input_reads == 0:
            warnings.warn(f"No reads counted in {reads}", UserWarning)

        # --- inputs that affect the signature ---
        self._hash_inputs = [reads, genome]

        # Conda environment and output prefix for CLI
        env        = cfg.run.conda_env
        out_prefix = f"{self.run_id}_flair"

        # --- parse flags from the config for this stage ---
        raw_flags = cfg.run.stages[0].flags  # first [[run.stages]] entry
        flag_parts: list[str] = []

        def _switch(k: str):
            # ensure flags begin with '--'
            flag_parts.append(k if k.startswith("--") else f"--{k}")

        if isinstance(raw_flags, list):
            # list of strings: use directly
            flag_parts.extend(raw_flags)

        elif isinstance(raw_flags, str):
            # space-separated string: split on whitespace
            flag_parts.extend(raw_flags.split())

        else:
            # TOML table: iterate attributes
            for key, val in vars(raw_flags).items():
                if val is False:
                    # user explicitly disabled this flag
                    continue
                if val in ("", None, True):
                    # boolean or switch flag (no value)
                    _switch(key)
                elif isinstance(val, (int, float)):
                    # numeric flag: add key and value
                    _switch(key)
                    flag_parts.append(str(val))
                else:
                    # file path flag: resolve relative to root/data_dir
                    p = Path(val)
                    if not p.is_absolute():
                        p = (Path(root) / data_dir / p).resolve()
                    _switch(key)
                    flag_parts.append(str(p))
                    # include file in signature if it exists, warn otherwise
                    if p.exists():
                        self._hash_inputs.append(p)
                    else:
                        warnings.warn(f"Flag file {p} does not exist; skipping it", UserWarning)

        # store flags for signature computation
        self._flags_components = flag_parts

        # --- capture `flair --version` once per run to include in signature ---
        if not hasattr(self, "_tool_version"):
            try:
                raw = subprocess.check_output(
                    ["conda", "run", "-n", env, "flair", "--version"],
                    text=True
                ).strip()
                # take last non-empty line
                self._tool_version = raw.splitlines()[-1] if raw else "flair-unknown"
            except subprocess.CalledProcessError:
                warnings.warn("Could not run `flair --version`; using 'flair-unknown'", UserWarning)
                self._tool_version = "flair-unknown"

        # --- debug printouts to stdout for visibility ---
        print("[DEBUG] align _hash_inputs:", self._hash_inputs)
        print("[DEBUG] align flags       :", " ".join(flag_parts))

        # --- warn if no extra flags were supplied ---
        if not flag_parts:
            warnings.warn("No extra flags configured for align stage; using defaults", UserWarning)

        # --- return the final command list for subprocess.call() ---
        return [
            "conda", "run", "-n", env,
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
