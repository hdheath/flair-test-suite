# src/flair_test_suite/stages/base.py
# -----------------------------------
# This abstract base class defines the common execution logic for all stages
# in the FLAIR test suite. Each stage (e.g., align, correct, slice) extends
# this class by implementing build_cmd() and expected_outputs().

from __future__ import annotations

import json         # to read/write marker files
import platform     # to record host info in markers
import shlex        # to split command strings
import subprocess   # to invoke external tools
import time         # to measure runtime
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List
import warnings

# core plumbing:
# PathBuilder: constructs the stage directory and filenames
# compute_signature, write_marker: create a unique signature and record metadata
# hash_many: compute a hash of input files for signature
# Reinstate: decide whether to skip, QC-only, or run
from ..lib.paths      import PathBuilder
from ..lib.signature  import compute_signature, write_marker
from ..lib.input_hash import hash_many
from ..lib.reinstate  import Reinstate

# QC lives in its own top-level package:
# QC_REGISTRY: maps stage names to collector functions
# qc_sidecar_path: constructs path to <stage>_qc.tsv
# load_marker: read existing marker JSON
from ..qc import QC_REGISTRY, qc_sidecar_path, load_marker


class StageBase(ABC):
    """
    Abstract base class for pipeline stages.

    Subclasses must implement:
      - build_cmd(): assemble command-line invocation
      - expected_outputs(): map output names to Path objects

    During build_cmd(), subclasses should populate:
      - self._hash_inputs: list of Paths/files affecting the signature
      - self._flags_components: list of flags for signature
    """
    # Unique name for registry lookups; override in each subclass
    name: str
    # Dependencies on other stages; override in each subclass if needed
    requires: tuple[str, ...] = ()
    # Which output key is considered the "primary" (e.g., bam, bed)
    primary_output_key: str = "bam"

    def __init__(
        self,
        cfg,
        run_id: str,
        work_dir: Path,
        upstreams: Dict[str, PathBuilder] | None = None,
    ):
        # Configuration namespace (parsed from TOML)
        self.cfg = cfg
        # Run identifier (used as a sub-directory)
        self.run_id = run_id
        # Root directory for all outputs
        self.work_dir = work_dir
        # PathBuilder objects for previously-run stages
        self.upstreams = upstreams or {}

    @abstractmethod
    def build_cmd(self) -> List[str] | str:
        """Build and return the command to run this stage"""
        ...

    @abstractmethod
    def expected_outputs(self) -> Dict[str, Path]:
        """Define expected output files for this stage"""
        ...

    def run(self):
        # 1) Let subclass assemble command and record inputs/flags
        raw_cmd = self.build_cmd()
        cmd = raw_cmd if isinstance(raw_cmd, list) else shlex.split(raw_cmd)

        # 2) Compute final hashes and flags string for signature
        self._input_hashes = hash_many(getattr(self, "_hash_inputs", []))
        self._flags_str = " ".join(getattr(self, "_flags_components", []))
        sig = self.signature  # compute or return existing

        # 3) Prepare the stage directory: work_dir/run_id/name/sig
        pb = PathBuilder(self.work_dir, self.run_id, self.name, sig)
        stage_dir = pb.stage_dir
        stage_dir.mkdir(parents=True, exist_ok=True)

        # Map expected outputs to their full paths
        outputs = {k: stage_dir / p for k, p in self.expected_outputs().items()}
        primary = outputs[self.primary_output_key]
        needs_qc = self.name in QC_REGISTRY

        # 4) Decide whether to skip, QC-only, or run the tool
        action = Reinstate.decide(
            stage_dir,
            primary,
            needs_qc=needs_qc,
            stage_name=self.name,
        )
        # Emit warnings based on decision
        if action == "run" and primary.exists():
            warnings.warn(
                f"Primary output exists before run for {self.name}; will overwrite",
                UserWarning
            )
        if action == "qc":
            warnings.warn(
                f"Primary exists but QC incomplete for {self.name}; regenerating QC",
                UserWarning
            )

        # --- skip case: tool already ran & QC exists ---
        if action == "skip":
            print(f"[SKIP] {self.name} complete (sig={sig})")
            marker_f = stage_dir / ".completed.json"
            try:
                full_meta = load_marker(marker_f)
                pb.metadata = full_meta.get("qc", {})
            except FileNotFoundError:
                pb.metadata = {}
            return pb

        # --- QC-only case: regenerate QC sidecar ---
        if action == "qc":
            print(f"[QC]   Regenerating QC for {self.name} (sig={sig})")
            qc_metrics = self._run_qc(stage_dir, primary, runtime=None)
            pb.metadata = qc_metrics
            marker_f = stage_dir / ".completed.json"
            old_meta = load_marker(marker_f)
            old_meta["qc"] = qc_metrics
            write_marker(pb, old_meta)
            return pb

        # --- run the external tool ---
        start = time.time()
        with open(stage_dir / "tool_stdout.log", "w") as out, \
             open(stage_dir / "tool_stderr.log", "w") as err:
            exit_code = subprocess.call(cmd, cwd=stage_dir, stdout=out, stderr=err)
            if exit_code != 0:
                warnings.warn(
                    f"External tool for stage {self.name} exited with code {exit_code}",
                    UserWarning
                )
        runtime = round(time.time() - start, 2)
        if not primary.exists():
            warnings.warn(
                f"Stage {self.name} completed but primary output missing: {primary}",
                UserWarning
            )

        # --- after run, collect QC if needed ---
        qc_metrics = self._run_qc(stage_dir, primary, runtime)
        pb.metadata = qc_metrics

        # --- write the completion marker with metadata ---
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

    def _run_qc(
        self,
        stage_dir: Path,
        primary: Path,
        runtime: float | None
    ) -> Dict:
        """
        Internal helper to invoke the registered QC collector function.
        Returns metrics dict or empty dict.
        """
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

    @property
    def signature(self) -> str:
        """
        Compute or return the stage signature based on tool, flags, and inputs.
        """
        if not hasattr(self, "_sig"):
            self._sig = compute_signature(
                self.tool_version,
                self.flags_str,
                self.input_hashes,
            )
        return self._sig

    @property
    def tool_version(self) -> str:
        """Default tool version; override in subclasses if needed."""
        return "flair-unknown"

    @property
    def flags_str(self) -> str:
        """Stringified flags used for signature calculation."""
        return getattr(self, "_flags_str", "")

    @property
    def input_hashes(self) -> List[str]:
        """List of file hashes used for signature calculation."""
        return getattr(self, "_input_hashes", [])
