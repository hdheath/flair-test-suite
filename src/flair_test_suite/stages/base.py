from __future__ import annotations

import logging
import platform
import shlex
import subprocess
import time
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from enum import Enum, auto
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from ..lib.paths import PathBuilder
from ..lib.signature import compute_signature, write_marker, qc_sidecar_path, load_marker
from ..lib.input_hash import hash_many
from ..lib.reinstate import Reinstate
from ..qc import QC_REGISTRY
from .stage_utils import get_stage_config, parse_cli_flags


class StageAction(Enum):
    RUN = auto()
    SKIP = auto()
    QC_ONLY = auto()


class StageBase(ABC):
    """
    Abstract base for pipeline stages.

    Subclasses must implement:
      • build_cmds()
      • expected_outputs()

    During build_cmd(s), subclasses should populate:
      • self._hash_inputs: list[Path|str] used for signature
      • self._flags_components: list[str] used for signature
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
        # memoized later:
        self._sig: Optional[str] = None
        self._flags_str: str = ""
        self._input_hashes: List[str] = []

    # ─────────────────────────────── Required API ───────────────────────────────

    @abstractmethod
    def build_cmds(self) -> List[List[str]] | List[str] | str:
        """Return command list(s) for the stage."""
        ...

    @abstractmethod
    def expected_outputs(self) -> Dict[str, Path]:
        ...

    # ─────────────────────────────── Public helpers ─────────────────────────────

    def resolve_stage_inputs(self, inputs: dict[str, object]) -> dict[str, object]:
        """
        Resolve stage inputs relative to cfg.run.data_dir and warn if missing.

        Accepts values that are a single path (str/Path) or a list/tuple of such
        values. Returns a dict mapping the same keys to resolved Path or
        list[Path].
        """
        from .stage_utils import resolve_path

        data_dir = Path(self.cfg.run.data_dir)
        resolved: dict[str, object] = {}
        for key, raw in inputs.items():
            # Normalize None
            if raw is None:
                resolved[key] = None
                continue

            # Handle iterables (list/tuple) of path-like values
            if isinstance(raw, (list, tuple)):
                paths = []
                for item in raw:
                    p = resolve_path(item, data_dir=data_dir)
                    paths.append(p)
                    if not p.exists():
                        logging.warning(f"[{self.name}] Input file missing: {key} -> {p}")
                resolved[key] = paths
            else:
                p = resolve_path(raw, data_dir=data_dir)
                resolved[key] = p
                if not p.exists():
                    logging.warning(f"[{self.name}] Input file missing: {key} -> {p}")
        return resolved

    def resolve_stage_flags(self, raw_flags=None) -> tuple[list[str], list[Path]]:
        """Turn a flags dict into CLI parts and collect extra file inputs for signature."""
        cfg = self.cfg
        data_dir = Path(cfg.run.data_dir)
        if raw_flags is None:
            stage_cfg = get_stage_config(cfg, self.name)
            raw_flags = getattr(stage_cfg, "flags", {}) or {}
        flag_parts, extra_inputs = parse_cli_flags(raw_flags, data_dir=data_dir)
        return flag_parts, extra_inputs

    # ─────────────────────────────── Main runner ────────────────────────────────

    def run(self) -> PathBuilder:
        """Main stage execution: build commands, run, QC, and write marker."""

        # 1) Build and normalize commands
        cmds = self._build_and_normalize_cmds()

        # 2) Prepare signature and stage directory
        self._prepare_signature()
        pb, outputs, primary, needs_qc = self._prepare_stage_dir()

        # 3) Handle no-work case
        if not cmds:
            return self._complete_empty_stage(pb)

        # 4) Decide action (run / skip / qc)
        action = self._decide_action(pb.stage_dir, primary, needs_qc)

        if action is StageAction.SKIP:
            return self._handle_skip(pb)
        if action is StageAction.QC_ONLY:
            return self._handle_qc_only(pb, primary)

        # 5) Execute commands
        runtime_sec, exit_code, started_ts = self._execute_commands(cmds, pb.stage_dir)

        # 6) Validate outputs
        self._validate_primary_output(primary)

        # 7) Run QC and finalize
        return self._finalize_stage(pb, cmds, runtime_sec, exit_code, started_ts, primary)

    # ────────── Helpers for run() orchestration ──────────

    def _prepare_stage_dir(self) -> Tuple[PathBuilder, Dict[str, Path], Path, bool]:
        pb = PathBuilder(self.work_dir, self.run_id, self.name, self.signature)
        stage_dir = pb.stage_dir
        stage_dir.mkdir(parents=True, exist_ok=True)

        outputs = {k: stage_dir / p for k, p in self.expected_outputs().items()}
        primary = outputs[self.primary_output_key]
        needs_qc = self.name in QC_REGISTRY
        return pb, outputs, primary, needs_qc

    def _decide_action(self, stage_dir: Path, primary: Path, needs_qc: bool) -> StageAction:
        decision = Reinstate.decide(
            stage_dir, primary, needs_qc=needs_qc, stage_name=self.name
        )
        if decision == "run":
            if primary.exists():
                logging.warning(
                    f"Primary output exists before run for {self.name}; will overwrite"
                )
            result = StageAction.RUN
        elif decision == "skip":
            result = StageAction.SKIP
        elif decision == "qc":
            logging.warning(
                f"Primary exists but QC incomplete for {self.name}; regenerating QC"
            )
            result = StageAction.QC_ONLY
        else:
            raise ValueError(f"Unknown Reinstate decision: {decision}")
        logging.info(f"[{self.name}] action: {result}")
        return result

    def _handle_skip(self, pb: PathBuilder) -> PathBuilder:
        logging.info(f"[SKIP] {self.name} complete (sig={self.signature})")
        try:
            full_meta = load_marker(pb.stage_dir)
            pb.metadata = full_meta if isinstance(full_meta, dict) else {}
        except FileNotFoundError:
            pb.metadata = {}
        return pb

    def _handle_qc_only(self, pb: PathBuilder, primary: Path) -> PathBuilder:
        logging.info(f"[QC] Regenerating QC for {self.name} (sig={self.signature})")
        qc_metrics = self._run_qc(pb.stage_dir, primary, runtime=None)
        pb.metadata = qc_metrics
        old_meta = load_marker(pb.stage_dir / ".completed.json")
        if not isinstance(old_meta, dict):
            old_meta = {}
        old_meta["qc"] = qc_metrics
        write_marker(pb, old_meta)
        return pb

    def _execute_commands(
        self, cmds: List[List[str]], stage_dir: Path
    ) -> Tuple[float, int, float]:
        start_ts = time.time()
        exit_code = self._run_all(cmds, log_path=stage_dir / "tool.log", cwd=stage_dir)
        runtime_sec = round(time.time() - start_ts, 2)
        logging.info(
            f"[{self.name}] commands finished in {runtime_sec}s (exit {exit_code})"
        )
        return runtime_sec, exit_code, start_ts

    def _validate_primary_output(self, primary: Path) -> None:
        if not primary.exists():
            logging.warning(f"Stage {self.name} completed but primary output missing: {primary}")
            raise RuntimeError(f"Stage {self.name} failed: missing primary output {primary}")

    def _finalize_stage(
        self,
        pb: PathBuilder,
        cmds: List[List[str]],
        runtime_sec: float,
        exit_code: int,
        started_ts: float,
        primary: Path,
    ) -> PathBuilder:
        qc_metrics = self._run_qc(pb.stage_dir, primary, runtime_sec)
        pb.metadata = qc_metrics

        meta = self._make_marker(
            cmds=[" ".join(c) for c in cmds],
            exit_code=exit_code,
            runtime_sec=runtime_sec,
            started_ts=started_ts,
            qc_metrics=qc_metrics,
        )
        write_marker(pb, meta)
        return pb

    # ─────────────────────────────── Internal helpers ───────────────────────────

    def _build_and_normalize_cmds(self) -> List[List[str]]:
        """Normalize build_cmds() output into a list of command lists."""
        raw_cmds = self.build_cmds()
        if raw_cmds is None:
            raw_cmds = []
        elif isinstance(raw_cmds, str):
            raw_cmds = [shlex.split(raw_cmds)]
        elif raw_cmds and isinstance(raw_cmds[0], (str, Path)):
            raw_cmds = [list(raw_cmds)]
        return [c if isinstance(c, list) else shlex.split(c) for c in raw_cmds]

    def _prepare_signature(self) -> None:
        """Memoize input hashes and flags string for signature computation."""
        self._input_hashes = hash_many(getattr(self, "_hash_inputs", []))
        self._flags_str = " ".join(getattr(self, "_flags_components", []))
        _ = self.signature  # trigger memoization

    def _complete_empty_stage(self, pb: PathBuilder) -> PathBuilder:
        """Write minimal QC + marker when there is nothing to run."""
        logging.info(f"[{self.name}] No commands to run; completing with minimal QC.")
        stage_dir = pb.stage_dir
        qc_summary = {"regions": 0, "qc_runtime_sec": 0.0, "all_skipped_empty": True}

        if self.name in QC_REGISTRY:
            from ..qc import write_metrics

            write_metrics(stage_dir, self.name, qc_summary)

        now = self._now_utc_iso()
        meta = {
            "stage": self.name,
            "signature": self.signature,
            "cmds": [],
            "exit_code": 0,
            "runtime_sec": 0.0,
            "started": now,
            "ended": now,
            "host": platform.node(),
            "qc": qc_summary,
        }
        write_marker(pb, meta)
        pb.metadata = qc_summary
        return pb

    def _run_all(self, cmds: List[List[str]], log_path: Path, cwd: Path) -> int:
        """Run all commands sequentially; raise on first non-zero exit."""
        exit_code = 0
        for cmd in cmds:
            exit_code = self.run_tool(cmd, log_path=log_path, cwd=cwd)
            if exit_code != 0:
                logging.warning(
                    f"External tool for stage {self.name} exited with code {exit_code}"
                )
                break
        return exit_code

    def _run_qc(self, stage_dir: Path, primary: Path, runtime: float | None) -> Dict:
        """Generic QC runner. Stages may override this (e.g., CorrectStage)."""
        logging.info(f"Executing _run_qc for stage: {self.name}")
        qc_func = QC_REGISTRY.get(self.name)
        if not qc_func or not primary.exists():
            return {}

        try:
            genome_fa = getattr(self, "_genome_fa_abs", None)
            kwargs = {
                "out_dir": stage_dir,
                "n_input_reads": getattr(self, "_n_input_reads", None),
                "runtime_sec": runtime,
            }
            if hasattr(self, "_read_count_method"):
                kwargs["read_count_method"] = getattr(self, "_read_count_method")
            if self.name == "align":
                metrics = qc_func(
                    primary,
                    genome_fa=genome_fa,
                    **kwargs,
                )
            else:
                metrics = qc_func(primary, **kwargs)
            logging.info(f"QC metrics for {self.name}: {metrics}")
            return metrics
        except Exception as e:
            logging.info(f"QC for '{self.name}' failed: {e}")
            return {}

    def _make_marker(
        self,
        cmds: List[str],
        exit_code: int,
        runtime_sec: float,
        started_ts: float,
        qc_metrics: Dict,
    ) -> Dict:
        """Create the marker dict consistently in one place."""
        meta = {
            "stage": self.name,
            "signature": self.signature,
            "cmds": cmds,
            "exit_code": exit_code,
            "runtime_sec": runtime_sec,
            "started": datetime.fromtimestamp(started_ts, timezone.utc).isoformat(),
            "ended": self._now_utc_iso(),
            "host": platform.node(),
            "qc": qc_metrics,
        }
        if hasattr(self, "_n_input_reads"):
            meta["n_input_reads"] = self._n_input_reads
        if hasattr(self, "_read_count_method"):
            meta["read_count_method"] = self._read_count_method
        return meta

    @staticmethod
    def _now_utc_iso() -> str:
        return datetime.now(timezone.utc).isoformat()

    # ───────────────────────────── Properties / I/O ─────────────────────────────

    @property
    def signature(self) -> str:
        if self._sig is None:
            self._sig = compute_signature(
                self.tool_version, self.flags_str, self.input_hashes
            )
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

    # ───────────────────────────── Subprocess wrapper ───────────────────────────

    def run_tool(self, cmd: list[str], log_path: Path = None, cwd: Path = None) -> int:
        """
        Run a single external command inside the configured conda env
        (unless the command already starts with 'conda').
        Logs to file, raises on non-zero exit.
        """
        env = self.cfg.run.conda_env
        full_cmd = ["conda", "run", "-n", env] + cmd if cmd and cmd[0] != "conda" else cmd
        log_path = log_path or Path("tool.log")
        logging.info(f"[{self.name}] Running: {' '.join(full_cmd)}")

        with open(log_path, "w") as logf:
            proc = subprocess.run(full_cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=cwd)

        if proc.returncode != 0:
            logging.error(
                f"[{self.name}] Command failed: {' '.join(full_cmd)} (exit {proc.returncode})"
            )
            raise RuntimeError(f"{self.name} failed with exit code {proc.returncode}")
        return proc.returncode

