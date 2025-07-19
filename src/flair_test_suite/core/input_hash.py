from __future__ import annotations
import hashlib, json
from pathlib import Path
from typing import Iterable, List, Union

__all__ = ["hash_path", "hash_scalar", "hash_many"]

def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def hash_path(p: Path) -> str:
    """sha256 of file *contents* (raises if file missing)."""
    return _sha256_bytes(p.read_bytes())

def hash_scalar(x: Union[str, int, float, bool]) -> str:
    """sha256 of a JSON‑encoded scalar (order‑invariant for bool/int/str)."""
    return _sha256_bytes(json.dumps(x, sort_keys=True).encode())

def hash_many(items: Iterable[Union[Path, str, int, float, bool]]) -> List[str]:
    """Convenience: hash a mixed iterable, returning **sorted‑unique** list."""
    out: set[str] = set()
    for it in items:
        if isinstance(it, Path):
            out.add(hash_path(it.resolve()))
        else:
            out.add(hash_scalar(it))
    return sorted(out)
