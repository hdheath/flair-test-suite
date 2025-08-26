# ────────────────────────────────────────────────────────────
# src/flair_test_suite/core/reinstate.py
# ------------------------------------------------------------------
from __future__ import annotations

import json
from pathlib import Path
from typing import Literal
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
               stage_name: str) -> Action:
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

        if needs_qc and stage_name in {"collapse", "transcriptome"}:
            ted_tsv = stage_dir / "qc" / "ted" / "TED.tsv"
            # Check for regionalized run by presence of region-tagged isoform beds
            is_regionalized = any(stage_dir.glob("*_*_*.isoforms.bed"))
            if is_regionalized:
                browser_dir = stage_dir / "qc" / "ted" / "transcriptome_browser"
                region_map = browser_dir / "region_map.json"
                if region_map.exists():
                    try:
                        mapping = json.loads(region_map.read_text())
                        browser_ok = all(Path(p).exists() for p in mapping.values())
                    except Exception:
                        browser_ok = False
                else:
                    browser_ok = False
                qc_ok = ted_tsv.exists() and browser_ok
                qc_desc = f"TED.tsv {'found' if ted_tsv.exists() else 'missing'}, browser plot {'found' if browser_ok else 'missing'}"
            else:
                qc_ok = ted_tsv.exists()
                qc_desc = f"TED.tsv {'found' if ted_tsv.exists() else 'missing'}"
        elif needs_qc:
            qc_tsv = qc_sidecar_path(stage_dir, stage_name)
            qc_ok = qc_tsv.exists()
            qc_desc = f"{qc_tsv} exists? {qc_ok}"
        else:
            qc_ok = True
            qc_desc = "n/a"

        print(f"[Reinstate] Stage: {stage_name}")
        print(f"  Marker: {marker_f} exists? {marker_ok}")
        print(f"  Primary: {primary} exists? {primary_ok}")
        print(f"  QC: {qc_desc}")
        print(f"  All files in {stage_dir}: {list(stage_dir.iterdir())}")

        if marker_ok and primary_ok and qc_ok:
            return "skip"

        if primary_ok and not qc_ok and needs_qc:
            return "qc"

        return "run"

