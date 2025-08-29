# ────────────────────────────────────────────────────────────
# src/flair_test_suite/core/reinstate.py
# ------------------------------------------------------------------
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Literal, Dict
from .signature import qc_sidecar_path, load_marker

# The three actions our state machine can return
Action = Literal["skip", "qc", "run"]


class Reinstate:
    """
    Stateless helper used by StageBase.run() to decide whether
    we can skip a stage, regenerate only its QC, or need to run it.
    """

    # ------------------------------------------------------------------
    @staticmethod
    def decide(stage_dir: Path,
               primary: Path,
               needs_qc: bool,
               stage_name: str,
               qc_expected: Dict[str, Path] | None = None) -> Action:
        """
        Parameters
        ----------
        stage_dir   : directory for this stage+signature
        primary     : primary output file (Path)
        needs_qc    : does this stage have a QC collector?
        stage_name  : plain name ("align", "correct", …)

        Returns
        -------
        "skip"  – marker, primary output *and* QC (if needed) exist
        "qc"    – primary output exists, QC missing ⇒ regenerate QC only
        "run"   – otherwise run the external tool
        """
        marker_f = stage_dir / ".completed.json"
        marker_ok = marker_f.exists()
        primary_ok = primary.exists()

        if needs_qc:
            # If the stage provided an explicit set of expected QC files, check all of them.
            if qc_expected:
                missing = [lbl for lbl, p in qc_expected.items() if not Path(p).exists()]
                qc_ok = len(missing) == 0
                if qc_ok:
                    qc_desc = f"all expected QC outputs present ({len(qc_expected)})"
                else:
                    qc_desc = "missing: " + ", ".join(missing)
            else:
                # Fall back to a single sidecar TSV
                qc_tsv = qc_sidecar_path(stage_dir, stage_name)
                qc_ok = qc_tsv.exists()
                qc_desc = f"{qc_tsv} exists? {qc_ok}"
        else:
            qc_ok = True
            qc_desc = "n/a"

        logger = logging.getLogger(__name__)
        logger.info(f"[Reinstate] Stage: {stage_name}")
        logger.info(f"  Marker: {marker_f} exists? {marker_ok}")
        logger.info(f"  Primary: {primary} exists? {primary_ok}")
        logger.info(f"  QC: {qc_desc}")
        try:
            files = [p.name for p in stage_dir.iterdir()]
        except Exception:
            files = []
        logger.info(f"  All files in {stage_dir}: {files}")

        # If the completion marker is missing but the primary exists and only QC
        # is missing, allow QC-only reinstatement when we have explicit QC expectations.
        if not marker_ok:
            if primary_ok and needs_qc and not qc_ok and qc_expected:
                return "qc"
            return "run"

        if marker_ok and primary_ok and qc_ok:
            return "skip"

        if primary_ok and not qc_ok and needs_qc:
            return "qc"

        return "run"
