from __future__ import annotations
import hashlib
from pathlib import Path 
import os


class PathBuilder:
    """Create stable output/layout paths."""

    def __init__(self, work_dir: Path, sample: str, stage: str, signature: str):
        self.base = Path(work_dir).expanduser().resolve()
        self.sample = sample
        self.stage = stage
        self.signature = signature

    # ---------- canonical locations ----------
    @property
    def stage_dir(self) -> Path:
        return self.base / self.sample / self.stage / self.signature

    def out(self, *parts: str | Path) -> Path:
        return self.stage_dir.joinpath(*map(str, parts))

    def marker(self) -> Path:
        return self.stage_dir / ".completed.json"

    # ---------- static helpers --------------
    @staticmethod
    def sha256(fp: Path, bs: int = 2**20) -> str:
        h = hashlib.sha256()
        with open(fp, "rb") as f:
            for chunk in iter(lambda: f.read(bs), b""):
                h.update(chunk)
        return h.hexdigest()

    @staticmethod
    def sha256_str(txt: str) -> str:
        import hashlib
        return hashlib.sha256(txt.encode()).hexdigest()

