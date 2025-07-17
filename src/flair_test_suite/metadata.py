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


def is_complete(pb: PathBuilder) -> bool:
    return pb.marker().is_file()
