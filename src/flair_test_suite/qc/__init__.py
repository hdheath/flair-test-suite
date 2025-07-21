"""
QC plug‑in registry  +  helpers used by StageBase.
KEEP ALL COLLECTOR MODULES (align_qc.py, correct_qc.py, …) IN THIS PACKAGE.
"""

from __future__ import annotations
import json
from pathlib import Path
from typing import Callable, Dict

# ------------------------------------------------------------- #
# Registry
# ------------------------------------------------------------- #
QC_REGISTRY: Dict[str, Callable] = {}

def register(stage_name: str):
    """Decorator: @register("align") auto‑registers the collector."""
    def _wrap(func: Callable):
        QC_REGISTRY[stage_name] = func
        return func
    return _wrap

# ------------------------------------------------------------- #
# Helper paths
# ------------------------------------------------------------- #
def qc_sidecar_path(stage_dir: Path, stage_name: str) -> Path:
    return stage_dir / f"{stage_name}_qc.tsv"

def load_marker(stage_dir: Path):
    m = stage_dir / ".completed.json"
    if not m.exists():
        return None
    try:
        return json.loads(m.read_text())
    except Exception:
        return None

# (Optional) tiny helper to write a “metric\tvalue” TSV in one line
def write_metrics(path: Path, metrics: Dict[str, object]):
    path.write_text(
        "metric\tvalue\n" + "\n".join(f"{k}\t{v}" for k, v in metrics.items())
    )
