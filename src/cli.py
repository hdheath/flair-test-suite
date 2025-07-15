## src/cli.py
import os
import sys
import logging
from .config_loader import load_config, Config


def main():
    if len(sys.argv) != 2:
        print("Usage: flair-test-suite <config.toml>", file=sys.stderr)
        sys.exit(1)

    config_path = sys.argv[1]
    raw_cfg      = load_config(config_path)
    cfg          = Config(raw_cfg)

    # Setup logging
    log_dir = 'logs'
    os.makedirs(log_dir, exist_ok=True)
    logging.basicConfig(
        filename=os.path.join(log_dir, f"run_{cfg.run_id}.log"),
        level=logging.ERROR,
        format="%(asctime)s %(levelname)s: %(message)s"
    )

    # Print summary and paths
    print(f"Run {cfg.run_id}: {cfg.summary}")
    dd = cfg.data_dir
    d  = cfg.data
    print(f"Dataset dir: {dd}")

    # Pipeline execution
    align_all(cfg)
    correct_all(cfg)
    slice_all_regions(cfg)
    collapse_all_regions(cfg)
    parse_all_regions(cfg)

    print("Completed FLAIR test suite.")
