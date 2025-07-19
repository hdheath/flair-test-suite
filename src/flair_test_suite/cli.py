# src/flair_test_suite/cli.py
import click, sys
from pathlib import Path
from .config_loader import load_config
from .stages import STAGE_REGISTRY
from .core import (              #  ← one dot, not two
    PathBuilder,
    topological_sort,
)

from typing import Dict

@click.command()
@click.argument("config", type=click.Path(exists=True, path_type=Path))
def main(config):
    cfg = load_config(config)
    sample = cfg.sample or cfg.run.sample          # no CLI override needed
    if not sample:
        raise click.UsageError("sample must be set in the TOML file")

    work_dir = Path(cfg.run.work_dir)
    stage_order = topological_sort(cfg.run.stages)

    upstreams: Dict[str, PathBuilder] = {}
    for st_cfg in stage_order:
        StageCls = STAGE_REGISTRY[st_cfg.name]
        pb = StageCls(cfg, sample, work_dir, upstreams).run()
        upstreams[st_cfg.name] = pb
        print(f"✓ Done {st_cfg.name} – outputs: {pb.stage_dir}")

if __name__ == "__main__":
    sys.exit(main())



