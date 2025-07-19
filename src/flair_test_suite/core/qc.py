"""
core.qc
=======

• Central registry where stages look up their QC collector.
• Decorator ``@register(stage_name)`` for collector modules.
• Helper functions for TSV side‑cars & marker editing.
"""
from __future__ import annotations
from pathlib import Path
import json
from typing import Callable, Dict, Any

# ------------------------------------------------------------------ #
# Registry + decorator
# ------------------------------------------------------------------ #
QC_REGISTRY: Dict[str, Callable[..., Dict[str, Any]]] = {}


def register(stage: str):
    """Decorator – add function to QC_REGISTRY under *stage* name."""
    def deco(fn):
        QC_REGISTRY[stage] = fn
        return fn
    return deco


# ------------------------------------------------------------------ #
# Helper paths / I/O
# ------------------------------------------------------------------ #
def qc_sidecar_path(stage_dir: Path, stage_name: str) -> Path:
    """Return <stage_name>_qc.tsv inside *stage_dir*."""
    return stage_dir / f"{stage_name}_qc.tsv"


def write_metrics(stage_dir: Path, stage_name: str, metrics: Dict[str, Any]):
    """Write keyed metrics as TSV and return written path."""
    tsv = qc_sidecar_path(stage_dir, stage_name)
    with tsv.open("w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in metrics.items():
            fh.write(f"{k}\t{v}\n")
    return tsv


def load_marker(stage_dir: Path) -> Dict[str, Any] | None:
    m = stage_dir / ".completed.json"
    if m.exists():
        return json.loads(m.read_text())
    return None

