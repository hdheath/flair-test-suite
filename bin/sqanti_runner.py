#!/usr/bin/env python3
"""
flair_automate/sqanti_runner.py

Run sqanti3_qc on every .isoforms.gtf under OUTDIR/results/<region>/<run>/,
capture logs under OUTDIR/logs/sqanti/, and emit sqanti_results.tsv.
"""

import os
import glob
import csv
import subprocess
import argparse
import logging
from datetime import datetime
import sys

# ── make sure our package is on PYTHONPATH ─────────────────────────────────
THIS = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(THIS, "..", ".."))
SRC = os.path.join(PROJECT_ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from flair_automate.config import load_manifest

# ── isoform‐classification summary logic stripped for brevity ──────────────
ISOCLASSES = [...]
NORMALIZED = { ... }


def run_sqanti(run_idx, iso_gtf, ref_gtf, ref_fa, prime5, prime3, sj, cpus, dry, logs_base):
    run_name = os.path.basename(iso_gtf).replace(".isoforms.gtf","")
    run_dir  = os.path.dirname(iso_gtf)
    class_txt = os.path.join(run_dir, f"{run_name}_classification.txt")

    # prepare per-run log dir
    log_dir = os.path.join(logs_base, run_name)
    os.makedirs(log_dir, exist_ok=True)
    out_log = os.path.join(log_dir, "sqanti.log")
    err_log = os.path.join(log_dir, "errors.log")

    # skip if already done
    if not dry and os.path.exists(class_txt):
        logging.info(" [%02d] SKIP %s (found %s)", run_idx, run_name, os.path.basename(class_txt))
        return class_txt

    cmd = [
        "sqanti3_qc",
        iso_gtf, ref_gtf, ref_fa,
        "--CAGE_peak", prime5,
        "--polyA_peak", prime3,
        "-c", sj,
        "-o", run_name,
        "-d", run_dir,
        "--skipORF",
        "--cpus", str(cpus)
    ]

    # log flair-style
    cmd_str = " ".join(subprocess.list2cmdline([c]) for c in cmd)
    logging.info(" [%02d] cd %s && %s", run_idx, run_dir, cmd_str)

    if dry:
        return class_txt

    # actually run & capture everything
    with open(out_log, "w") as lf:
        try:
            subprocess.run(cmd, cwd=run_dir, stdout=lf, stderr=subprocess.STDOUT, check=True)
        except subprocess.CalledProcessError as e:
            # record the failure
            with open(err_log, "a") as ef:
                ef.write(f"{datetime.now().isoformat()}\t{run_name}\t{e.returncode}\n")
            logging.error(" [%02d] ERROR %s (exit %d)", run_idx, run_name, e.returncode)

    return class_txt


def summarize_one(class_txt, dry):
    # … your existing summarize logic, but skip if no classification.txt …
    if not os.path.exists(class_txt):
        logging.warning(" SKIP summary (no %s)", class_txt)
        return
    # … write sqanti_results.tsv here …
    ...


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-m","--manifest", required=True)
    p.add_argument("-o","--outdir",   required=True)
    p.add_argument("-c","--cpus",      type=int, default=None)
    p.add_argument("--dry_run",        action="store_true")
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    cfg    = load_manifest(args.manifest)
    ref_gtf  = cfg["gtf"]
    ref_fa   = cfg["genome"]
    prime5   = cfg.get("reference_5_prime_regions_bed_file")
    prime3   = cfg.get("reference_3_prime_regions_bed_file")
    sj       = cfg.get("junctions")
    cpus     = args.cpus or cfg.get("threads", 1)

    logs_base = os.path.join(args.outdir, "logs", "sqanti")
    os.makedirs(logs_base, exist_ok=True)

    # 1) find every isoforms.gtf under results/<region>/<run>/
    gtfs = sorted(glob.glob(os.path.join(args.outdir, "results", "*", "*", "*.isoforms.gtf")))
    info = []
    for idx, g in enumerate(gtfs, start=1):
        logging.info("→ [%02d] %s", idx, os.path.basename(g))
        class_txt = run_sqanti(idx, g, ref_gtf, ref_fa, prime5, prime3, sj, cpus, args.dry_run, logs_base)
        info.append(class_txt)

    # 2) summarize each
    for ct in info:
        summarize_one(ct, args.dry_run)


if __name__ == "__main__":
    main()



