from __future__ import annotations

import platform, shlex, subprocess, time
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

from ..core import PathBuilder, QC_REGISTRY
from ..core.input_hash import hash_many
from ..core.reinstate import Reinstate, Action
from ..core.signature import compute_signature, write_marker


class StageBase(ABC):
    """Common scaffold inherited by every pipeline stage."""

    # subclasses override ------------------------------------------------
    name: str
    requires: tuple[str, ...] = ()
    primary_output_key: str = "bam"

    # -------------------------------------------------------------------
    def __init__(
        self,
        cfg,
        sample: str,
        work_dir: Path,
        upstreams: Dict[str, PathBuilder] | None = None,
    ):
        self.cfg        = cfg
        self.sample     = sample
        self.work_dir   = work_dir
        self.upstreams  = upstreams or {}

        # child stages simply append → self._hash_inputs / _flags_components
        self._hash_inputs: List[Path] = []
        self._flags_components: List[str] = []

    # ---------- subclass hooks -----------------------------------------
    @abstractmethod
    def build_cmd(self) -> List[str] | str: ...

    @abstractmethod
    def expected_outputs(self) -> Dict[str, Path]: ...

    # -------------------------------------------------------------------
    def run(self):
        # 1) let subclass build the command & populate inputs / flags
        raw_cmd = self.build_cmd()
        cmd = raw_cmd if isinstance(raw_cmd, list) else shlex.split(raw_cmd)

        # 2) canonicalise signature pieces
        self._input_hashes = hash_many(self._hash_inputs)
        self._flags_str    = " ".join(self._flags_components)

        sig = self.signature
        pb  = PathBuilder(self.work_dir, self.sample, self.name, sig)
        stage_dir = pb.stage_dir
        stage_dir.mkdir(parents=True, exist_ok=True)

        # 3) ask Reinstate what to do
        outputs  = {k: stage_dir / p for k, p in self.expected_outputs().items()}
        primary  = outputs[self.primary_output_key]
        action   = Reinstate.decide(stage_dir, primary, needs_qc=self.name in QC_REGISTRY)

        # --- SKIP -------------------------------------------------------
        if action is Action.SKIP:
            print(f"[SKIP] {self.name} complete (sig={sig})")
            return pb

        # --- QC‑only ----------------------------------------------------
        if action is Action.QC_ONLY:
            print(f"[QC]   Regenerating QC for {self.name} (sig={sig})")
            self._run_qc(stage_dir, primary, runtime=None)
            return pb

        # --- RUN tool ---------------------------------------------------
        start = time.time()
        with open(stage_dir / "tool_stdout.log", "w") as out, \
             open(stage_dir / "tool_stderr.log", "w") as err:
            exit_code = subprocess.call(cmd, cwd=stage_dir, stdout=out, stderr=err)
        runtime = round(time.time() - start, 2)

        qc_metrics = self._run_qc(stage_dir, primary, runtime)

        # marker ---------------------------------------------------------
        meta = {
            "stage":        self.name,
            "signature":    sig,
            "cmd":          " ".join(cmd),
            "exit_code":    exit_code,
            "runtime_sec":  runtime,
            "started":      datetime.fromtimestamp(start, timezone.utc).isoformat(),
            "ended":        datetime.now(timezone.utc).isoformat(),
            "host":         platform.node(),
            "qc":           qc_metrics,
        }
        if hasattr(self, "_n_input_reads"):
            meta["n_input_reads"] = self._n_input_reads
        write_marker(pb, meta)
        return pb

    # ---------- QC helper ----------------------------------------------
    def _run_qc(self, stage_dir: Path, primary: Path, runtime: float | None) -> Dict:
        qc_fn = QC_REGISTRY.get(self.name)
        if not qc_fn or not primary.exists():
            return {}
        try:
            return qc_fn(
                primary,
                out_dir       = stage_dir,
                n_input_reads = getattr(self, "_n_input_reads", None),
                runtime_sec   = runtime,
            )
        except Exception as e:
            print(f"[WARN] QC for '{self.name}' failed: {e}")
            return {}

    # ---------- signature props ----------------------------------------
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
    def tool_version(self) -> str:
        return "flair-unknown"

    @property
    def flags_str(self) -> str:        # used by compute_signature
        return self._flags_str

    @property
    def input_hashes(self) -> List[str]:
        return self._input_hashes
