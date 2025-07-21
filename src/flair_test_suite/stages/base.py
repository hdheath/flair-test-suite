# ── src/flair_test_suite/stages/base.py ───────────────────────────────
from __future__ import annotations

import json, platform, shlex, subprocess, time
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

# core plumbing
from ..core.paths        import PathBuilder
from ..core.signature    import compute_signature, write_marker
from ..core.input_hash   import hash_many
from ..core.reinstate    import Reinstate            # ← central restart logic

# QC lives in its own top‑level package
from ..qc                import QC_REGISTRY, qc_sidecar_path


# ──────────────────────────────────────────────────────────────────────
class StageBase(ABC):
    """
    Shared execution scaffold for all pipeline stages.

    Sub‑classes must implement:
        · build_cmd()          -> list[str] | str
        · expected_outputs()   -> dict[str, Path]

    During build_cmd() they should populate:
        · self._hash_inputs      (iterable[Path | str])
        · self._flags_components (iterable[str])
    """
    # to override in concrete stage
    name: str
    requires: tuple[str, ...] = ()
    primary_output_key: str = "bam"

    # ──────────────────────────────────────────────────────────────────
    def __init__(
        self,
        cfg,
        sample: str,
        work_dir: Path,
        upstreams: Dict[str, PathBuilder] | None = None,
    ):
        self.cfg       = cfg
        self.sample    = sample
        self.work_dir  = work_dir
        self.upstreams = upstreams or {}

    # ===== sub‑class API =============================================
    @abstractmethod
    def build_cmd(self) -> List[str] | str: ...
    @abstractmethod
    def expected_outputs(self) -> Dict[str, Path]: ...

    # ===== driver ====================================================
    def run(self):
        # 1. Let subclass assemble command (and fill hash/flag pieces)
        raw_cmd = self.build_cmd()
        cmd     = raw_cmd if isinstance(raw_cmd, list) else shlex.split(raw_cmd)

        # 2. Finalise signature inputs & flags
        self._input_hashes = hash_many(getattr(self, "_hash_inputs", []))
        self._flags_str    = " ".join(getattr(self, "_flags_components", []))
        sig = self.signature

        # 3. Stage directory
        pb         = PathBuilder(self.work_dir, self.sample, self.name, sig)
        stage_dir  = pb.stage_dir
        stage_dir.mkdir(parents=True, exist_ok=True)

        outputs   = {k: stage_dir / p for k, p in self.expected_outputs().items()}
        primary   = outputs[self.primary_output_key]
        needs_qc  = self.name in QC_REGISTRY

        # 4. Ask Reinstate what to do
        action = Reinstate.decide(
            stage_dir,
            primary,
            needs_qc=(self.name in QC_REGISTRY),
            stage_name=self.name,             # ← pass the right name
        )

        if action == "skip":
            print(f"[SKIP] {self.name} complete (sig={sig})")
            return pb

        if action == "qc":
            print(f"[QC]   Regenerating QC for {self.name} (sig={sig})")
            self._run_qc(stage_dir, primary, runtime=None)
            return pb

        # 5. Run external tool
        start = time.time()
        with open(stage_dir / "tool_stdout.log", "w") as out, \
             open(stage_dir / "tool_stderr.log", "w") as err:
            exit_code = subprocess.call(cmd, cwd=stage_dir, stdout=out, stderr=err)
        runtime = round(time.time() - start, 2)

        # 6. QC after run
        qc_metrics = self._run_qc(stage_dir, primary, runtime)

        # 7. Marker
        meta = {
            "stage":      self.name,
            "signature":  sig,
            "cmd":        " ".join(cmd),
            "exit_code":  exit_code,
            "runtime_sec": runtime,
            "started":    datetime.fromtimestamp(start, timezone.utc).isoformat(),
            "ended":      datetime.now(timezone.utc).isoformat(),
            "host":       platform.node(),
            "qc":         qc_metrics,
        }
        if hasattr(self, "_n_input_reads"):
            meta["n_input_reads"] = self._n_input_reads
        write_marker(pb, meta)
        return pb

    # ── QC helper ─────────────────────────────────────────────────────
    def _run_qc(self,
                stage_dir: Path,
                primary:   Path,
                runtime:   float | None) -> Dict:
        qc_func = QC_REGISTRY.get(self.name)
        if not qc_func or not primary.exists():
            return {}

        try:
            return qc_func(
                primary,
                out_dir=stage_dir,
                n_input_reads=getattr(self, "_n_input_reads", None),
                runtime_sec=runtime,
            )
        except Exception as e:
            print(f"[WARN] QC for '{self.name}' failed: {e}")
            return {}

    # ── signature properties ─────────────────────────────────────────
    @property
    def signature(self) -> str:
        if not hasattr(self, "_sig"):
            self._sig = compute_signature(
                self.tool_version,
                self.flags_str,
                self.input_hashes,
            )
        return self._sig

    @property
    def tool_version(self) -> str:     # overridden in subclasses
        return "flair-unknown"

    @property
    def flags_str(self) -> str:
        return getattr(self, "_flags_str", "")

    @property
    def input_hashes(self) -> List[str]:
        return getattr(self, "_input_hashes", [])
