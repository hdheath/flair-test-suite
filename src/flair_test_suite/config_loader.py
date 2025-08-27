import tomllib
from pathlib import Path
from .config_schema import Config
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

    for st in stage_paths:
        with open(st, "rb") as fh:
            frag = tomllib.load(fh)
        frag_run = frag.get("run", {})
        frag_stages = frag_run.get("stages")
        if not frag_stages:
            raise ValueError(f"Stage config {st} has no [run.stages] entries")
        stages.extend(frag_stages)

    return Config.parse_obj(data)

