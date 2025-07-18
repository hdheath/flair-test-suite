import click, sys
from pathlib import Path
from .config_loader import load_config
from .stages import STAGE_REGISTRY


@click.command()
@click.argument("config", type=click.Path(exists=True, path_type=Path))
@click.option("--sample", help="Sample name (overrides config)")
def main(config: Path, sample: str | None):
    cfg = load_config(config)

    # precedence: CLI flag > top‑level key > [run].sample
    sample = (
        sample
        or getattr(cfg, "sample", None)
        or getattr(cfg.run, "sample", None)
    )
    if not sample:
        raise click.UsageError("No sample specified in CLI or config")

    work_dir = Path(cfg.run.work_dir)

    for stage_cfg in cfg.run.stages:
        Stage = STAGE_REGISTRY[stage_cfg.name]
        pb = Stage(cfg, sample, work_dir).run()
        if pb:
            print(f"✓ Done {stage_cfg.name} – outputs: {pb.stage_dir}")


if __name__ == "__main__":
    sys.exit(main())


