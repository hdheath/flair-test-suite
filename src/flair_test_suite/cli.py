#!/usr/bin/env python3
import argparse
import sys

from flair_test_suite.config import load_config
from flair_test_suite.printer import print_config_summary

def main():
    parser = argparse.ArgumentParser(
        prog="flair-test-suite",
        description="Run the FLAIR Test Suite: load & print your TOML config"
    )

    # Positional, required argument:
    parser.add_argument(
        "config",
        metavar="CONFIG_TOML",
        help="Path to your TOML config file (e.g. config/config_template.toml)"
    )

    args = parser.parse_args()

    try:
        cfg = load_config(args.config)
    except Exception as e:
        print(f"[error] failed to load config '{args.config}': {e}", file=sys.stderr)
        sys.exit(1)

    print_config_summary(cfg)

if __name__ == "__main__":
    main()
