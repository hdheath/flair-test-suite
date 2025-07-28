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
import logging

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
from .lib.logging import setup_run_logging  # <-- centralized logging

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

    # Set up temporary logging to terminal only
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.StreamHandler()]
    )

    any_ran = False
    any_failed = False
    all_configs_skipped = False 
    seen_run_ids = set()  # Track run_ids for this invocation

    for cfg_path in inputs:
        if not cfg_path.exists():
            warnings.warn(f"Config file not found: {cfg_path}", UserWarning)
            continue

        logging.info(f"Loading config: {cfg_path}")

        try:
            cfg = load_config(cfg_path)
        except Exception as e:
            warnings.warn(f"Failed to load {cfg_path}: {e}", UserWarning)
            continue

        run_id = getattr(cfg, "run_id", None) or getattr(cfg.run, "run_id", None)
        work_dir = Path(cfg.run.work_dir)
        log_path = work_dir / run_id / "run_summary.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)

        # Decide log mode: overwrite if first time in this invocation, append otherwise
        mode = "w" if run_id not in seen_run_ids else "a"
        seen_run_ids.add(run_id)

        # Centralized logging setup for this run_id
        setup_run_logging(log_path, mode)

        # Topologically sort stages for this config
        stage_order = topological_sort(cfg.run.stages)

        upstreams: Dict[str, PathBuilder] = {}
        logging.info(f"Executing run '{run_id}' with {len(stage_order)} stage(s)")
        all_skipped = True

        for st_cfg in stage_order:
            StageCls = STAGE_REGISTRY[st_cfg.name]
            stage_instance = StageCls(cfg, run_id, work_dir, upstreams)
            stage_instance.build_cmd()
            try:
                pb = stage_instance.run()
            except Exception as e:
                logging.error(f"Stage {st_cfg.name} failed: {e}")
                any_failed = True
                break
            else:
                upstreams[st_cfg.name] = pb
                logging.info(f"✓ Done {st_cfg.name} – outputs: {pb.stage_dir}")
                if getattr(stage_instance, "action", None) != "skip":
                    all_skipped = False
                    any_ran = True
        if not all_skipped:
            all_configs_skipped = False

    # After the loop:
    if seen_run_ids and all_configs_skipped:
        logging.info("All configurations are up to date; nothing to do.")
        sys.exit(0)
    elif not seen_run_ids:
        raise click.UsageError("No valid configuration files were executed. Check the input paths.")
    elif any_failed:
        logging.error("Pipeline failed: at least one stage did not complete successfully.")
        sys.exit(1)
    elif any_ran:
        logging.info("All stages completed successfully.")

        if not all_skipped:
            all_configs_skipped = False  # <-- Add this after each config

if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
