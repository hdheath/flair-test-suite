from __future__ import annotations
import sys, click
from pathlib import Path

from .config_loader import load_config
from .stages import STAGE_REGISTRY


@click.command()
@click.argument("config", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--sample", default="SAMPLE", help="Sample identifier")
def main(config: Path, sample: str):
    cfg = load_config(config)
    work_dir = Path(cfg.run.work_dir)

    stage_cfg = cfg.run.stages[0]          # we’re only running align for now
    Stage     = STAGE_REGISTRY[stage_cfg.name]

    pb = Stage(cfg, sample, work_dir).run()
    print(f"✓ Done {stage_cfg.name} – outputs: {pb.stage_dir}")


if __name__ == "__main__":
    sys.exit(main())

