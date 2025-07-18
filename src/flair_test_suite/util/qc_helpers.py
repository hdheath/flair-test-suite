from __future__ import annotations
from pathlib import Path
import json

def qc_sidecar_path(stage_dir: Path, stage_name: str) -> Path:
    """Return <stage_name>_qc.tsv in the stage directory."""
    return stage_dir / f"{stage_name}_qc.tsv"

def load_marker(stage_dir: Path) -> dict | None:
    """Read .completed.json if present; else return None."""
    m = stage_dir / ".completed.json"
    if m.exists():
        return json.loads(m.read_text())
    return None
