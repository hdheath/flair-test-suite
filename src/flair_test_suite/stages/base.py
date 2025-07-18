from __future__ import annotations
import subprocess, shlex, time, platform
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from pathlib import Path
import json
from ..paths import PathBuilder
from ..metadata import compute_signature, write_marker, is_complete
from ..qc import QC_REGISTRY             # <─ new central registry
from ..util.qc_helpers import qc_sidecar_path, load_marker

class StageBase(ABC):
    name: str
    requires: tuple[str, ...] = ()
    primary_output_key: str = "bam"

    def __init__(
        self,
        cfg,
        sample: str,
        work_dir: Path,
        upstreams: dict[str, "PathBuilder"] | None = None,  # ← new
    ):
        self.cfg        = cfg
        self.sample     = sample
        self.work_dir   = work_dir
        self.upstreams  = upstreams or {}


    # ------------------------------------------------------------------ #
    # Methods subclasses
    # ------------------------------------------------------------------ #
    @abstractmethod
    def build_cmd(self) -> list[str] | str:
        """Return the command to execute (list preferred, raw string allowed)."""

    @abstractmethod
    def expected_outputs(self) -> dict[str, Path]:
        """Return key → Path mapping of main output files."""

    # ------------------------------------------------------------------ #
    # Main execution
    # ------------------------------------------------------------------ #
    def run(self):
        # prerequisite check omitted for brevity …

        cmd = self.build_cmd() if isinstance(self.build_cmd(), list) else shlex.split(self.build_cmd())

        pb = PathBuilder(self.work_dir, self.sample, self.name, self.signature)
        stage_dir = pb.stage_dir
        stage_dir.mkdir(parents=True, exist_ok=True)

        # ---------- fast‑path: outputs already present --------------------
        expected = [stage_dir / p for p in self.expected_outputs().values()]
        qc_tsv   = qc_sidecar_path(stage_dir, self.name)
        marker   = load_marker(stage_dir)

        outputs_ok = all(p.exists() for p in expected)
        qc_ok      = qc_tsv.exists() and marker and marker.get("qc")

        if outputs_ok and qc_ok:
            print(f"[SKIP] {self.name} already complete (sig={self.signature})")
            return pb

        # ---------- (B) outputs ok, QC missing ---------------------------
        if outputs_ok and not qc_ok and self.name in QC_REGISTRY:
            print(f"[QC]   Regenerating QC for {self.name} (sig={self.signature})")
            qc_metrics = QC_REGISTRY[self.name](
                expected[0],          # primary output
                out_dir=stage_dir,
                n_input_reads=marker.get("n_input_reads") if marker else None,
                runtime_sec=marker.get("runtime_sec") if marker else None,
            )
            if marker:
                marker["qc"] = qc_metrics
                (stage_dir / ".completed.json").write_text(
                    json.dumps(marker, indent=2)
                )
            return pb

        # ---------- (C) need to run the tool -----------------------------
        start = time.time()
        with open(stage_dir / "tool_stdout.log", "w") as out, \
            open(stage_dir / "tool_stderr.log", "w") as err:
            exit_code = subprocess.call(cmd, cwd=stage_dir, stdout=out, stderr=err)

        runtime = time.time() - start


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
