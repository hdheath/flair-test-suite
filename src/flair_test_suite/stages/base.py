from __future__ import annotations
import subprocess, shlex, time, platform
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from pathlib import Path

from ..paths import PathBuilder
from ..metadata import compute_signature, write_marker, is_complete
from ..qc import QC_REGISTRY             # <─ new central registry


class StageBase(ABC):
    """Common execution scaffold for all pipeline stages (align, correct, …)."""

    name: str  # must be overridden in each concrete subclass
    primary_output_key: str = "bam"   # <─ NEW default for QC look‑up

    # ------------------------------------------------------------------ #
    def __init__(self, cfg, sample: str, work_dir: Path):
        self.cfg = cfg
        self.sample = sample
        self.work_dir = work_dir

    # ------------------------------------------------------------------ #
    # Methods subclasses MUST implement
    # ------------------------------------------------------------------ #
    @abstractmethod
    def build_cmd(self) -> list[str] | str:
        """Return the command to execute (list preferred, raw string allowed)."""

    @abstractmethod
    def expected_outputs(self) -> dict[str, Path]:
        """Return key → Path mapping of main output files."""

    # ------------------------------------------------------------------ #
    # Main execution engine
    # ------------------------------------------------------------------ #
    def run(self):
        cmd = self.build_cmd()
        if isinstance(cmd, str):                       # allow raw strings
            cmd = shlex.split(cmd)

        sig = self.signature
        pb = PathBuilder(self.work_dir, self.sample, self.name, sig)

        # ---------------------------------------------------------- skip logic
        if (pb.stage_dir / ".completed.json").exists():
            exp_files = self.expected_outputs().values()
            missing = [fp for fp in exp_files if not (pb.stage_dir / fp).exists()]

            # QC side‑car required if a collector is registered
            qc_sidecar = None
            if self.name in QC_REGISTRY:
                qc_sidecar = pb.stage_dir / f"{self.name}_qc.tsv"
                if not qc_sidecar.exists():
                    missing.append(qc_sidecar)

            if not missing:
                print(f"[SKIP] {self.name} already complete (sig={sig})")
                return pb       # ← early return only if *everything* is present
            else:
                print(f"[WARN] Re‑running {self.name}: missing {len(missing)} files")


        pb.stage_dir.mkdir(parents=True, exist_ok=True)

        # ---------- execute ----------
        start = time.time()
        with open(pb.stage_dir / "tool_stdout.log", "w") as out, \
            open(pb.stage_dir / "tool_stderr.log", "w") as err:
            exit_code = subprocess.call(
                cmd,
                cwd=pb.stage_dir,
                stdout=out,
                stderr=err,
            )

        runtime   = round(time.time() - start, 2)

        # ---------- QC collection --------------------------------------------
        qc_metrics = {}
        qc_func = QC_REGISTRY.get(self.name)
        if qc_func:
            key = getattr(self, "primary_output_key", "bam")
            primary = self.expected_outputs().get(key)
            if primary and not primary.is_absolute():
                primary = pb.stage_dir / primary

            if primary.exists():
                try:
                    qc_metrics = qc_func(
                        primary,                         # ← positional, no keyword
                        out_dir=pb.stage_dir,
                        n_input_reads=getattr(self, "_n_input_reads", None),
                        runtime_sec=runtime,
                    )
                except Exception as e:
                    print(f"[WARN] QC for stage '{self.name}' failed: {e}")
            else:
                print(f"[WARN] QC skipped: {primary} does not exist")

        # ---------- write completion marker ----------------------------------
        meta = {
            "stage": self.name,
            "signature": sig,
            "cmd": " ".join(map(str, cmd)),
            "exit_code": exit_code,
            "runtime_sec": round(runtime, 2),
            "started": datetime.fromtimestamp(start, timezone.utc).isoformat(),
            "ended":   datetime.now(timezone.utc).isoformat(),
            "host": platform.node(),
            "qc": qc_metrics,
        }
        if hasattr(self, "_n_input_reads"):
            meta["n_input_reads"] = self._n_input_reads
        write_marker(pb, meta)
        return pb


    # ------------------------------------------------------------------ #
    # Signature helpers
    # ------------------------------------------------------------------ #
    @property
    def signature(self) -> str:
        if not hasattr(self, "_sig"):
            self._sig = compute_signature(
                self.tool_version,
                self.flags_str,
                self.input_hashes,
            )
        return self._sig

    # --- sensible fall-backs; subclasses override as needed ------------
    @property
    def tool_version(self) -> str:
        return "flair-unknown"

    @property
    def flags_str(self) -> str:
        return ""

    @property
    def input_hashes(self) -> list[str]:
        return []
