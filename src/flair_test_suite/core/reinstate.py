"""
Central logic for deciding whether a stage should be
 • skipped                      (all good)
 • re‑QC only                   (outputs exist, QC missing)
 • fully re‑run                 (outputs missing / stale)

The object is intentionally tiny so it can be reused by every stage.
"""

from __future__ import annotations
from enum import Enum, auto
from pathlib import Path
import json

from . import qc_sidecar_path


class Action(Enum):
    SKIP     = auto()
    QC_ONLY  = auto()
    RUN      = auto()


class Reinstate:
    """Stateless helper – all decisions funnel through one place."""

    @staticmethod
    def _marker_ok(marker_f: Path) -> bool:
        if not marker_f.exists():
            return False
        try:
            _ = json.loads(marker_f.read_text())
            return True
        except Exception:
            return False

    # ------------------------------------------------------------------
    @staticmethod
    def decide(
        stage_dir: Path,
        primary_out: Path,
        needs_qc: bool = False,
    ) -> Action:
        """
        Decide what to do for this stage instance.

        Parameters
        ----------
        stage_dir    : Path to the stage signature folder
        primary_out  : Primary output file (Path)
        needs_qc     : If True, the stage registers a QC collector

        Returns
        -------
        Action.SKIP     – everything present (incl. QC if required)
        Action.QC_ONLY  – outputs present but QC missing
        Action.RUN      – tool must be executed
        """
        marker_f = stage_dir / ".completed.json"

        if not primary_out.exists() or not Reinstate._marker_ok(marker_f):
            return Action.RUN

        if not needs_qc:
            return Action.SKIP

        qc_f = qc_sidecar_path(stage_dir, stage_dir.name)
        try:
            qc_block = json.loads(marker_f.read_text()).get("qc")
        except Exception:
            qc_block = None

        if qc_f.exists() and qc_block:
            return Action.SKIP

        return Action.QC_ONLY

