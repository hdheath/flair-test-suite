#!/usr/bin/env python3
"""
src/flair_test_suite/cli.py

Command-line interface for the FLAIR test suite:
- Load one or more TOML configs
- Set up logging per run
- Print config summary
- Execute alignment stage
- Run align pytest 
"""
import argparse
import sys
import logging
from pathlib import Path
import subprocess

from flair_test_suite.config_loader import load_config
from flair_test_suite.print_config import print_config_summary
from flair_test_suite.align import run_align


def setup_logging(run_id: str):
    """
    Prepare logging: create logs directory, set up a file and console handler.
    """
    logs_dir = Path("logs")
    logs_dir.mkdir(exist_ok=True)
    logfile = logs_dir / f"{run_id}.log"
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    # remove any existing handlers
    for handler in list(root.handlers):
        root.removeHandler(handler)
    fmt = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    # file handler
    fh = logging.FileHandler(logfile, mode="w")
    fh.setFormatter(fmt)
    root.addHandler(fh)
    # console (stdout) handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)
    root.addHandler(ch)


def main():
    parser = argparse.ArgumentParser(
        prog="flair-test-suite",
        description="FLAIR test suite CLI: run pipeline and auto-tests"
    )
    parser.add_argument(
        "configs", metavar="CONFIG_TOML", nargs='+',
        help="Paths to one or more TOML config files to run"
    )
    args = parser.parse_args()

    for cfg_path in args.configs:
        try:
            cfg = load_config(cfg_path)
        except Exception as e:
            print(
                f"[error] Failed to load config '{cfg_path}': {e}",
                file=sys.stderr
            )
            continue

        setup_logging(run_id=str(cfg.run_id))
        logging.info(f"=== Starting run {cfg.run_id} using {cfg_path} ===")
        print_config_summary(cfg)

        # Alignment stage
        if hasattr(cfg, 'run') and hasattr(cfg.run.steps, 'align'):
            try:
                run_align(cfg)
            except Exception:
                logging.exception(
                    f"Error during alignment for run {cfg.run_id}"
                )

            # Automatically run align tests after align stage
            logging.info(f"Running align tests for sample {cfg.sample_name}")
            test_ret = subprocess.call(["pytest", "tests/align", "-q"])
            if test_ret != 0:
                logging.error(
                    f"Align tests failed for sample {cfg.sample_name}; skipping this run"
                )
                continue

        logging.info(f"=== Finished run {cfg.run_id} ===\n")


if __name__ == "__main__":
    main()

