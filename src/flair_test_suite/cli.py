# src/flair_test_suite/cli.py
# --------------------------
# Modified CLI to accept a *TSV file* containing relative paths (one per line)
# to individual TOML configuration files. Each config is loaded and executed
# independently in the order they appear in the TSV.
#
# The TSV is expected to be a single column. Lines beginning with '#' or blank
# lines are ignored. Relative paths are resolved against the TSV file's parent
# directory. If a listed config file does not exist, a warning is printed and
# execution continues with the remaining configs.

import click   # CLI framework
import sys     # for sys.exit()
import warnings
from pathlib import Path  # for path manipulations

# load configuration from a TOML file
from .config_loader import load_config

# registry mapping stage names to their implementing classes
from .stages import STAGE_REGISTRY

# utility functions and classes for path building and DAG sorting
from .core import (
    PathBuilder,
    topological_sort,
)

from typing import Dict


def _iter_config_paths(tsv_path: Path):
    """Yield absolute Paths to config files listed in the TSV.

    Rules:
      * Skip blank lines and lines starting with '#'.
      * Use the first column before a tab (single-column TSV).
      * Resolve relative paths against the TSV's parent directory.
    """
    for raw in tsv_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        rel = line.split("\t")[0].strip()  # first column only
        cfg_path = Path(rel)
        if not cfg_path.is_absolute():
            cfg_path = (tsv_path.parent / cfg_path).resolve()
        yield cfg_path


@click.command()
@click.argument("config_list", type=click.Path(exists=True, path_type=Path))
def main(config_list: Path):
    """
    Main entrypoint for the CLI.

    CONFIG_LIST: Path to a single-column TSV file listing relative paths to
    TOML configuration files (one per line). Each config is executed in order.
    """
    any_ran = False  # track if at least one config succeeded

    for cfg_path in _iter_config_paths(config_list):
        if not cfg_path.exists():
            warnings.warn(f"Config file not found: {cfg_path}", UserWarning)
            continue

        print(f"[INFO] Loading config: {cfg_path}")
        try:
            cfg = load_config(cfg_path)
        except Exception as e:  # defensive: misformatted TOML, etc.
            warnings.warn(f"Failed to load {cfg_path}: {e}", UserWarning)
            continue

        # Determine sample name: prefer cfg.sample if set, else cfg.run.sample
        sample = getattr(cfg, "sample", None) or getattr(cfg.run, "sample", None)
        if not sample:
            warnings.warn(
                f"Skipping {cfg_path} – sample must be set in the TOML file.",
                UserWarning,
            )
            continue

        work_dir = Path(cfg.run.work_dir)
        # Topologically sort stages for this config
        stage_order = topological_sort(cfg.run.stages)

        upstreams: Dict[str, PathBuilder] = {}
        print(f"[INFO] Executing sample '{sample}' with {len(stage_order)} stage(s)")
        for st_cfg in stage_order:
            StageCls = STAGE_REGISTRY[st_cfg.name]
            pb = StageCls(cfg, sample, work_dir, upstreams).run()
            upstreams[st_cfg.name] = pb
            print(f"✓ Done {st_cfg.name} – outputs: {pb.stage_dir}")
        any_ran = True

    if not any_ran:
        raise click.UsageError(
            "No valid configuration files were executed. Check the TSV paths."
        )


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
