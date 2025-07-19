from __future__ import annotations
import json, hashlib, base64
from datetime import datetime, timezone
import platform
from .paths import PathBuilder    # instead of ..paths


def compute_signature(tool_version: str, flags: str, hashes: list[str]) -> str:
    uniq = sorted(set(hashes))               # stable order, de‑dup
    payload = f"{tool_version}|{flags}|{'|'.join(uniq)}"
    return hashlib.sha1(payload.encode()).hexdigest()[:12]



def write_marker(pb: PathBuilder, meta: dict):
    pb.stage_dir.mkdir(parents=True, exist_ok=True)
    with open(pb.marker(), "w") as fh:
        json.dump(meta, fh, indent=2)


def qc_sidecar_path(stage_dir: Path, stage_name: str) -> Path:
    return stage_dir / f"{stage_name}_qc.tsv"

def load_marker(stage_dir: Path) -> dict | None:
    m = stage_dir / ".completed.json"
    if m.exists():
        import json
        return json.loads(m.read_text())
    return None

def is_complete(stage_dir: Path, outputs: list[Path], needs_qc: bool) -> bool:
    """Return True if marker + outputs (+ optional QC) are present."""
    marker = stage_dir / ".completed.json"
    if not marker.exists():
        return False
    if not all(p.exists() for p in outputs):
        return False
    if needs_qc:
        qc_tsv = stage_dir / f"{stage_dir.name}_qc.tsv"
        if not qc_tsv.exists():
            return False
        if not json.loads(marker.read_text()).get("qc"):
            return False
    return True

