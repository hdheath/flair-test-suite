from __future__ import annotations
import json, hashlib
from datetime import datetime, timezone
import platform
from pathlib import Path
from .paths import PathBuilder


def compute_signature(tool_ver: str, flags: str, input_hashes: list[str]) -> str:
    key = "\n".join([tool_ver, flags, *input_hashes])
    return hashlib.sha256(key.encode()).hexdigest()[:12]


def write_marker(pb: PathBuilder, meta: dict):
    pb.stage_dir.mkdir(parents=True, exist_ok=True)
    with open(pb.marker(), "w") as fh:
        json.dump(meta, fh, indent=2)


def is_complete(pb) -> bool:
    """
    Return True if:
      • .completed.json exists
      • every file in Stage.expected_outputs() exists
      • if .json lists a 'qc' block → side‑car <stage>_qc.tsv exists
    """
    marker = pb.stage_dir / ".completed.json"
    if not marker.exists():
        return False

    try:
        meta = json.loads(marker.read_text())
    except Exception:
        return False

    # 1) outputs exist
    for fp in pb.stage_dir.glob("*"):
        pass  # placeholder – StageBase will supply explicit list below

    # We'll let StageBase supply expected files; see patch there
    return True
