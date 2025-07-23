# src/flair_test_suite/cli.py
# --------------------------
# Command-line interface for FLAIR test suite.
# Accepts either:
# 1) Path to a single TOML config: runs that config
# 2) Path to a single-column TSV listing multiple TOMLs

import click   # CLI framework
import sys     # for sys.exit()
import warnings
from pathlib import Path  # for path manipulations

# load configuration from a TOML file
from .config_loader import load_config

# registry mapping stage names to their implementing classes
from .stages import STAGE_REGISTRY

# utility functions and classes for path building and DAG sorting
from .lib import (
    PathBuilder,
    topological_sort,
)

from typing import Dict, Iterator


def _iter_config_paths(tsv_path: Path) -> Iterator[Path]:
    """Yield absolute Paths to config files listed in the TSV.

    Rules:
      * Skip blank lines and lines starting with '#'.
      * Use the first column before a tab (single-column TSV).
      * Resolve relative paths against the TSV file's parent directory.
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
@click.argument("config_input", type=click.Path(exists=True, path_type=Path))
def main(config_input: Path):
    """
    Main entrypoint for the CLI.

    CONFIG_INPUT may be either:
      - a path to a single TOML configuration file, or
      - a path to a single-column TSV file listing TOML config paths.

    If a TSV is provided, each listed TOML is loaded and executed in order.
    """
    # Determine if input is a list of configs or a single TOML
    if config_input.suffix.lower() in (".tsv", ".txt"):
        inputs = _iter_config_paths(config_input)
    else:
        inputs = iter([config_input])

    any_ran = False  # track if at least one config succeeded

    for cfg_path in inputs:
        if not cfg_path.exists():
            warnings.warn(f"Config file not found: {cfg_path}", UserWarning)
            continue

        print(f"[INFO] Loading config: {cfg_path}")
        try:
            cfg = load_config(cfg_path)
        except Exception as e:  # misformatted TOML, etc.
            warnings.warn(f"Failed to load {cfg_path}: {e}", UserWarning)
            continue

        # Determine run identifier: prefer cfg.run_id if set, else cfg.run.run_id
        run_id = getattr(cfg, "run_id", None) or getattr(cfg.run, "run_id", None)
        if not run_id:
            warnings.warn(
                f"Skipping {cfg_path} – run_id must be set in the TOML file.",
                UserWarning,
            )
            continue

        work_dir = Path(cfg.run.work_dir)
        # Topologically sort stages for this config
        stage_order = topological_sort(cfg.run.stages)

        upstreams: Dict[str, PathBuilder] = {}
        print(f"[INFO] Executing run '{run_id}' with {len(stage_order)} stage(s)")
        for st_cfg in stage_order:
            StageCls = STAGE_REGISTRY[st_cfg.name]
            pb = StageCls(cfg, run_id, work_dir, upstreams).run()
            upstreams[st_cfg.name] = pb
            print(f"✓ Done {st_cfg.name} – outputs: {pb.stage_dir}")
        any_ran = True

    if not any_ran:
        raise click.UsageError(
            "No valid configuration files were executed. Check the input paths."
        )


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())


