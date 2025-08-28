import tomllib
from pathlib import Path
from .config_schema import Config
from typing import Dict, Any


def _deep_merge(dst: Dict[str, Any], src: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge ``src`` into ``dst`` (in place) and return ``dst``.

    Only dict->dict merges are recursive; other values overwrite.
    """
    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            _deep_merge(dst[k], v)
        else:
            dst[k] = v
    return dst
from typing import Iterable

def load_config(path: str | Path) -> Config:
    with open(path, "rb") as fh:
        data = tomllib.load(fh)
    return Config.parse_obj(data)


def load_config_fragments(
    base_path: str | Path, stage_paths: Iterable[str | Path]
) -> Config:
    """Load a base config and append stages from ``stage_paths``."""
    with open(base_path, "rb") as fh:
        data = tomllib.load(fh)

    run = data.setdefault("run", {})
    stages = run.setdefault("stages", [])
    qc = data.setdefault("qc", {})

    for st in stage_paths:
        with open(st, "rb") as fh:
            frag = tomllib.load(fh)
        frag_run = frag.get("run", {})
        frag_stages = frag_run.get("stages")
        if not frag_stages:
            raise ValueError(f"Stage config {st} has no [run.stages] entries")
        stages.extend(frag_stages)

        # Merge top-level QC blocks from fragments so users can colocate
        # [qc.<stage>.*] with the stage template.
        if frag_qc := frag.get("qc"):
            _deep_merge(qc, frag_qc)

    return Config.parse_obj(data)
