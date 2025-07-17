#!/usr/bin/env python3
"""
src/flair_test_suite/cli.py

Command‑line interface for the FLAIR test suite:
- Load one or more TOML configs
- Set up logging per run
- Print config summary
- Execute each stage listed under [[run.stages]]
- Run pytest for each stage, if tests exist, passing the active config path
"""
import argparse
import sys
import logging
import subprocess
import os
from pathlib import Path
from importlib import import_module

from flair_test_suite.config_loader import load_config
from flair_test_suite.print_config import print_config_summary


def setup_logging(run_id: str):
    """Prepare logging: create logs directory, set up a file and console handler."""
    logs_dir = Path("logs")
    logs_dir.mkdir(exist_ok=True)
    logfile = logs_dir / f"{run_id}.log"

    root = logging.getLogger()
    root.setLevel(logging.INFO)
    for handler in list(root.handlers):
        root.removeHandler(handler)

    fmt = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")

    fh = logging.FileHandler(logfile, mode="w")
    fh.setFormatter(fmt)
    root.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)
    root.addHandler(ch)


def main():
    parser = argparse.ArgumentParser(
        prog="flair-test-suite",
        description="FLAIR test suite CLI: run pipeline and auto-tests",
    )
    parser.add_argument(
        "configs", metavar="CONFIG_TOML", nargs="+", help="Path(s) to TOML config(s)"
    )
    args = parser.parse_args()

    for cfg_path in args.configs:
        cfg_path_abs = str(Path(cfg_path).resolve())

        try:
            cfg = load_config(cfg_path_abs)
        except Exception as e:
            print(f"[error] Failed to load config '{cfg_path}': {e}", file=sys.stderr)
            continue

        setup_logging(run_id=str(cfg.run_id))
        logging.info(f"=== Starting run {cfg.run_id} using {cfg_path} ===")
        print_config_summary(cfg)

        # Carry the active run_id from one stage to the next
        active_run_id = str(cfg.run_id)
        for stage_cfg in cfg.run.stages:
            name = stage_cfg.name
            logging.info(f"--- Running stage: {name} (run_id={active_run_id}) ---")

            try:
                mod = import_module(f"flair_test_suite.stages.{name}")
                active_run_id = mod.run(stage_cfg.flags, cfg, active_run_id)
            except Exception:
                logging.exception(
                    f"Error during stage '{name}' for run {active_run_id}")
                break

            # Run pytest for this stage if tests/<stage> exists
            test_dir = Path("tests") / name
            if test_dir.exists():
                logging.info(f"Running pytest for stage '{name}'")
                env = os.environ.copy()
                env["FLAIR_CFG_PATH"] = cfg_path_abs  # pass config to tests
                ret = subprocess.call(["pytest", str(test_dir), "-q"], env=env)
                if ret != 0:
                    logging.error(
                        f"Tests failed for stage '{name}'; aborting run")
                    break

        logging.info(
            f"=== Finished run {cfg.run_id} (final run_id={active_run_id}) ===\n"
        )


if __name__ == "__main__":
    main()
