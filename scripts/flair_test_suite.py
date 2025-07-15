#!/usr/bin/env python3
"""
scripts/flair_test_suite.py: top-level entry point for per-region FLAIR pipeline.
Loads Python manifest, then runs FLAIR‐align, FLAIR‐correct, slicing,
FLAIR-collapse, and optionally SQANTI3 QC + plotting.
"""
import os
import sys
import argparse
import logging
import subprocess
import shlex
from importlib import util

# ─── make sure our package is importable ─────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT       = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
SRC        = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

# ─── imports from our library ────────────────────────────────────────────────
from flair_automate.slicer                     import slice_all_regions
from flair_automate.flair_collapse_automation           import collapse_all_regions
from flair_automate.parse_region_metrics_from_gtf import parse_all_regions
from flair_automate.flair_align_automation    import align_all
from flair_automate.flair_correct_automation  import correct_all


def get_args():
    p = argparse.ArgumentParser(
        description="Driver: run FLAIR align, correct, slice, collapse, (opt) SQANTI"
    )
    p.add_argument(
        '-m', '--manifest_py',
        default=os.path.join(ROOT, 'config', 'manifest.py'),
        help="Path to Python manifest module (default: config/manifest.py)"
    )
    p.add_argument(
        '-o', '--outdir',
        default=os.path.join(ROOT, 'outputs'),
        help="Base output directory (default: outputs/)"
    )
    p.add_argument(
        '--dry_run',
        action='store_true',
        help="Print commands without executing"
    )
    p.add_argument(
        '--run_sqanti',
        action='store_true',
        help="After collapse, run SQANTI3 QC and plotting"
    )
    return p.parse_args()


def main():
    args = get_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s"
    )

    # ── 1) Load pure-Python manifest ────────────────────────────────────────
    spec = util.spec_from_file_location("manifest", args.manifest_py)
    manifest = util.module_from_spec(spec)
    spec.loader.exec_module(manifest)
    cfg = manifest.CONFIG

    # ── 2) FLAIR‐ALIGN ───────────────────────────────────────────────────────
    logging.info("=== Running FLAIR align ===")
    align_all(cfg, args.outdir, dry_run=args.dry_run)

    # ── 3) FLAIR‐CORRECT ─────────────────────────────────────────────────────
    logging.info("=== Running FLAIR correct ===")
    correct_all(cfg, args.outdir, dry_run=args.dry_run)

    # ── 4) Slice every region ─────────────────────────────────────────────
    logging.info("=== Slicing all regions ===")
    slice_all_regions(cfg, args.outdir, dry_run=args.dry_run)

    # ── 5) Compute region GTF metrics ─────────────────────────────────────
    logging.info("=== Parsing region GTF metrics ===")
    parse_all_regions(
        os.path.join(args.outdir, "regions"),
        dry_run=args.dry_run
    )

    # ── 6) Run FLAIR collapse grid ────────────────────────────────────────
    logging.info("=== Running FLAIR collapse grid ===")
    collapse_all_regions(cfg, args.outdir, dry=args.dry_run)

'''
    # ── 7) Optional SQANTI3 QC + plotting ────────────────────────────────
    if args.run_sqanti:
        # a) SQANTI runner
        sq_script = os.path.join(SRC, "flair_automate", "sqanti_runner.py")
        sq_cmd = [
            sys.executable, sq_script,
            "--manifest_py", args.manifest_py,
            "--outdir",      args.outdir
        ]
        if args.dry_run:
            sq_cmd.append("--dry_run")
            logging.info("[DRY] %s", " ".join(shlex.quote(x) for x in sq_cmd))
        else:
            logging.info("[RUN] %s", " ".join(shlex.quote(x) for x in sq_cmd))
            subprocess.run(sq_cmd, check=True)

        # b) SQANTI plotting
        plot_script = os.path.join(SRC, "flair_automate", "sqanti_plot.py")
        plot_cmd = [sys.executable, plot_script, "--outdir", args.outdir]
        if args.dry_run:
            logging.info("[DRY] %s", " ".join(shlex.quote(x) for x in plot_cmd))
        else:
            logging.info("[RUN] %s", " ".join(shlex.quote(x) for x in plot_cmd))
            subprocess.run(plot_cmd, check=True)
'''

if __name__ == '__main__':
    main()

