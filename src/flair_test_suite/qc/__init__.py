# src/flair_test_suite/qc/__init__.py
# ----------------------------------
# QC plug‑in registry and helper functions for StageBase.
# All individual QC collectors (e.g., align_qc.py, correct_qc.py) live in this package.
# Ensure every collector module is imported so it can register itself.

from __future__ import annotations
from importlib import import_module  # dynamic import of plugin modules
from pathlib import Path  # for filesystem paths
import csv    # to write TSV metrics files
import json   # to read existing marker JSON files
import time   # optional timing utilities, not currently used directly

# ---------------------------------------------------------------------
# 1) Public API -------------------------------------------------------

# Global registry mapping stage names ("align", "correct", etc.) to their QC functions
QC_REGISTRY: dict[str, callable] = {}


def register(stage_name: str):
    """
    Decorator to register a QC collector function for a given stage.

    Usage:
      @register("align")
      def collect(...):
          return {...metrics...}
    """
    def _wrap(func):
        QC_REGISTRY[stage_name] = func
        return func
    return _wrap


def qc_sidecar_path(stage_dir: Path, stage_name: str) -> Path:
    """
    Construct the path to the <stage_name>_qc.tsv side‑car file inside stage_dir.
    """
    return stage_dir / f"{stage_name}_qc.tsv"


def write_metrics(stage_dir: Path, stage_name: str, metrics: dict):
    """
    Write the QC metrics dict to a TSV file under stage_dir.
    Filename is determined via qc_sidecar_path().

    TSV format:
      metric<tab>value
      key1<tab>val1
      key2<tab>val2
      ...
    """
    out_f = qc_sidecar_path(stage_dir, stage_name)
    with open(out_f, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        for k, v in metrics.items():
            w.writerow([k, v])


def load_marker(marker_f: Path) -> dict:
    """
    Read and parse the JSON completion marker file at marker_f.
    Returns the full metadata dict, including any "qc" block.
    """
    with open(marker_f) as fh:
        return json.load(fh)

# ---------------------------------------------------------------------
# 2) Auto‑import plug‑ins ---------------------------------------------
# Dynamically load each collector module so its @register decorator runs
for _name in ("align_qc", "correct_qc"):
    import_module(f"{__name__}.{_name}")

# Placeholder: add new QC collector modules here as needed
# e.g., import_module(f"{__name__}.collapse_qc")  
