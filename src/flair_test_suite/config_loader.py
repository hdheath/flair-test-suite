# src/flair_test_suite/config_loader.py
import tomli
from types import SimpleNamespace
from pathlib import Path


def _to_ns(obj):
    """Recursively convert dicts to SimpleNamespace, lists-of-dicts to list‑of‑NS."""
    if isinstance(obj, dict):
        return SimpleNamespace(**{k: _to_ns(v) for k, v in obj.items()})
    if isinstance(obj, list):
        return [_to_ns(v) for v in obj]
    return obj


def load_config(path: str | Path):
    with open(path, "rb") as fh:
        data = tomli.load(fh)
    return _to_ns(data)

