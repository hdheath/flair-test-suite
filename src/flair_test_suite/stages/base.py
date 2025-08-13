from __future__ import annotations

import platform
import shlex
import subprocess
import time
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List
import logging

from ..lib.paths      import PathBuilder
from ..lib.signature  import compute_signature, write_marker, qc_sidecar_path, load_marker
from ..lib.input_hash import hash_many
from ..lib.reinstate  import Reinstate
from ..qc import QC_REGISTRY
from .stage_utils import get_stage_config, parse_cli_flags


class StageBase(ABC):
    """
    Abstract base class for pipeline stages.
    Subclasses must implement:
      - build_cmds() (preferred) OR build_cmd() (legacy)
      - expected_outputs()
    During build_cmd(s), populate:
      - self._hash_inputs: list[Path/str] used for signature
      - self._flags_components: list[str] of flags used for signature
    """
    name: str
    requires: tuple[str, ...] = ()
    primary_output_key: str = "bam"

    def __init__(
        self,
        cfg,
        run_id: str,
        work_dir: Path,
        upstreams: Dict[str, PathBuilder] | None = None,
    ):
        self.cfg = cfg
        self.run_id = run_id
        self.work_dir = work_dir
        self.upstreams = upstreams or {}

    @abstractmethod
    def build_cmd(self) -> List[str] | str:
        """Legacy single-command builder. Prefer build_cmds()."""
        ...

    @abstractmethod
    def build_cmds(self) -> List[List[str]] | None:
        """Return a list of commands, or None to use build_cmd()."""
        return None

    def resolve_stage_inputs(self, inputs: dict[str, str | Path]) -> dict[str, Path]:
        data_dir = Path(self.cfg.run.data_dir)
        from .stage_utils import resolve_path
        resolved = {}
        for key, raw in inputs.items():
            p = resolve_path(raw, data_dir=data_dir)
            resolved[key] = p
            if not p.exists():
                logging.warning(f"[{self.name}] Input file missing: {key} -> {p}")
        return resolved

    def resolve_stage_flags(self, raw_flags=None) -> tuple[list[str], list[Path]]:
        cfg = self.cfg
        data_dir = Path(cfg.run.data_dir)
        if raw_flags is None:
            stage_cfg = get_stage_config(cfg, self.name)
            raw_flags = getattr(stage_cfg, "flags", {}) or {}
        flag_parts, extra_inputs = parse_cli_flags(raw_flags, data_dir=data_dir)
        return flag_parts, extra_inputs

    @abstractmethod
    def expected_outputs(self) -> Dict[str, Path]:
        ...

    def run(self):
        # Build commands (prefer build_cmds)
        raw_cmds = self.build_cmds()
        if raw_cmds is None:
            # legacy fallback to build_cmd()
            built = self.build_cmd()
            raw_cmds = [built] if isinstance(built, (list, str)) else []

        cmds = [cmd if isinstance(cmd, list) else shlex.split(cmd) for cmd in raw_cmds]

        # Compute signature inputs/flags (subclass must have set these)
        self._input_hashes = hash_many(getattr(self, "_hash_inputs", []))
        self._flags_str = " ".join(getattr(self, "_flags_components", []))
        sig = self.signature

        # Prepare stage directory
        pb = PathBuilder(self.work_dir, self.run_id, self.name, sig)
        stage_dir = pb.stage_dir
        stage_dir.mkdir(parents=True, exist_ok=True)

        # No work to do: write minimal QC and marker, return cleanly
        if not cmds:
            logging.info(f"[{self.name}] No commands to run; completing with minimal QC.")
            qc_summary = {"regions": 0, "qc_runtime_sec": 0.0, "all_skipped_empty": True}
            # If a collector is registered, give it a chance to write a stage-level sidecar name
            if self.name in QC_REGISTRY:
                from ..qc import write_metrics
                write_metrics(stage_dir, self.name, qc_summary)
            meta = {
                "stage": self.name,
                "signature": sig,
                "cmds": [],
                "exit_code": 0,
                "runtime_sec": 0.0,
                "started": datetime.now(timezone.utc).isoformat(),
                "ended":   datetime.now(timezone.utc).isoformat(),
                "host": platform.node(),
                "qc": qc_summary,
            }
            write_marker(pb, meta)
            return pb

        # Map expected outputs to absolute paths and find primary
        outputs = {k: stage_dir / p for k, p in self.expected_outputs().items()}
        primary = outputs[self.primary_output_key]
        needs_qc = self.name in QC_REGISTRY

        # Decide run/skip/qc
        action = Reinstate.decide(stage_dir, primary, needs_qc=needs_qc, stage_name=self.name)
        if action == "run" and primary.exists():
            logging.warning(f"Primary output exists before run for {self.name}; will overwrite")
        if action == "qc":
            logging.warning(f"Primary exists but QC incomplete for {self.name}; regenerating QC")

        if action == "skip":
            logging.info(f"[SKIP] {self.name} complete (sig={sig})")
            try:
                full_meta = load_marker(stage_dir)
                pb.metadata = full_meta if isinstance(full_meta, dict) else {}
            except FileNotFoundError:
                pb.metadata = {}
            return pb

        if action == "qc":
            logging.info(f"[QC] Regenerating QC for {self.name} (sig={sig})")
            qc_metrics = self._run_qc(stage_dir, primary, runtime=None)
            pb.metadata = qc_metrics
            # Update marker
            old_meta = load_marker(stage_dir / ".completed.json")
            if not isinstance(old_meta, dict):
                old_meta = {}
            old_meta["qc"] = qc_metrics
            write_marker(pb, old_meta)
            return pb

        # --- run tool(s) ---
        start = time.time()
        exit_code = 0
        for cmd in cmds:
            exit_code = self.run_tool(cmd, log_path=stage_dir / "tool.log", cwd=stage_dir)
            if exit_code != 0:
                logging.warning(f"External tool for stage {self.name} exited with code {exit_code}")
        runtime = round(time.time() - start, 2)

        if not primary.exists():
            logging.warning(f"Stage {self.name} completed but primary output missing: {primary}")
            raise RuntimeError(f"Stage {self.name} failed: missing primary output {primary}")

        # --- QC ---
        qc_metrics = self._run_qc(stage_dir, primary, runtime)
        pb.metadata = qc_metrics

        # --- marker ---
        meta = {
            "stage":      self.name,
            "signature":  sig,
            "cmds":       [" ".join(cmd) for cmd in cmds],
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

    def _run_qc(self, stage_dir: Path, primary: Path, runtime: float | None) -> Dict:
        """
        Generic QC runner. Stages may override this (e.g., CorrectStage).
        """
        logging.info(f"Executing _run_qc for stage: {self.name}")
        qc_func = QC_REGISTRY.get(self.name)
        if not qc_func:
            logging.warning(f"No QC function registered for stage: {self.name}")
            return {}
        if not primary.exists():
            logging.warning(f"Primary output does not exist for stage: {self.name}")
            return {}

        try:
            genome_fa = getattr(self, "_genome_fa_abs", None)
            if self.name == "align":
                return qc_func(
                    primary,
                    out_dir=stage_dir,
                    n_input_reads=getattr(self, "_n_input_reads", None),
                    genome_fa=genome_fa,
                    runtime_sec=runtime,
                )
            # Default pass-through for other stages with simple collectors
            return qc_func(
                primary,
                out_dir=stage_dir,
                n_input_reads=getattr(self, "_n_input_reads", None),
                runtime_sec=runtime,
            )
        except Exception as e:
            logging.info(f"QC for '{self.name}' failed: {e}")
            return {}

    @property
    def signature(self) -> str:
        if not hasattr(self, "_sig"):
            self._sig = compute_signature(self.tool_version, self.flags_str, self.input_hashes)
        return self._sig

    @property
    def tool_version(self) -> str:
        return "flair-unknown"

    @property
    def flags_str(self) -> str:
        return getattr(self, "_flags_str", "")

    @property
    def input_hashes(self) -> List[str]:
        return getattr(self, "_input_hashes", [])

    def run_tool(self, cmd: list[str], log_path: Path = None, cwd: Path = None) -> int:
        env = self.cfg.run.conda_env
        full_cmd = ["conda", "run", "-n", env] + cmd if cmd[0] != "conda" else cmd
        log_path = log_path or Path("tool.log")
        logging.info(f"[{self.name}] Running: {' '.join(full_cmd)}")
        with open(log_path, "w") as logf:
            proc = subprocess.run(full_cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=cwd)
        if proc.returncode != 0:
            logging.error(
                f"[{self.name}] Command failed: {' '.join(full_cmd)} "
                f"(exit {proc.returncode})"
            )
            raise RuntimeError(
                f"{self.name} failed with exit code {proc.returncode}"
            )
        return proc.returncode
