#!/usr/bin/env python3
"""
flair_automate/flair_collapse_automation.py

Run FLAIR collapse for each region, sample, and correction run,
using exactly the sets of flags you specified in your manifest
(no combinatorial expansion). Defaults to 4 threads if not set.
"""
import os
import glob
import shlex
import subprocess
import logging
from datetime import datetime

__all__ = ["flair_cmd", "run_jobs", "collapse_all_regions"]


def flair_cmd(rcfg):
    """
    Build the core 'flair collapse' command for one region (no -o, no flag args).
    rcfg must contain:
      - genome (str)
      - query_reg (path to .bed)
      - reads_fa_reg (path to .fasta)
      - gtf_reg   (path to .gtf)
      - threads   (int)
      - junctions_reg (optional path to .tab or None)
    """
    cmd = [
        "flair", "collapse",
        "-g", rcfg["genome"],
        "-q", rcfg["query_reg"],
        "-r", rcfg["reads_fa_reg"],
        "-f", rcfg["gtf_reg"],
        "-t", str(rcfg.get("threads"))
    ]
    return cmd


def run_jobs(flag_sets, base_cmd, results_base, logs_base, dry=False):
    """
    For each dict in flag_sets, run:
      cd results_base/<run_id> && /usr/bin/time -v flair collapse ... [flags] -o <run_id>
    and log to logs_base/<run_id>/collapse.time.log,
    skipping any run where that log already exists.
    """
    os.makedirs(results_base, exist_ok=True)
    os.makedirs(logs_base, exist_ok=True)
    err_log = os.path.join(logs_base, "errors.log")
    with open(err_log, "w") as ef:
        ef.write(f"# {datetime.now().isoformat()}\n")

    for i, flags in enumerate(flag_sets, start=1):
        # build stable run_id from the flags dict
        parts, args = [], []
        for flag, val in flags.items():
            if isinstance(val, bool):
                if val:
                    args.append(flag)
                    parts.append(flag.lstrip("-"))
            else:
                args += [flag, str(val)]
                parts.append(f"{flag.lstrip('-')}-{val}")
        run_id = "_".join(parts) or "default"

        res_dir = os.path.join(results_base, run_id)
        log_dir = os.path.join(logs_base, run_id)
        log_file = os.path.join(log_dir, "collapse.time.log")

        # skip if already done
        if not dry and os.path.exists(log_file):
            logging.info("    [SKIP] %s (already completed)", run_id)
            continue

        os.makedirs(res_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)

        # build full command
        full_cmd = base_cmd + args + ["-o", run_id]
        time_cmd = ["/usr/bin/time", "-v"] + full_cmd
        cmd_str  = " ".join(shlex.quote(x) for x in time_cmd)

        if dry:
            logging.info("    [DRY] cd %s && %s", res_dir, cmd_str)
        else:
            logging.info("    [%02d] cd %s && %s", i, res_dir, cmd_str)
            try:
                with open(log_file, "w") as lf:
                    subprocess.run(
                        time_cmd,
                        cwd=res_dir,
                        stdout=lf, stderr=subprocess.STDOUT,
                        check=True
                    )
            except subprocess.CalledProcessError as e:
                logging.error("    [ERROR] %s failed (exit %d)", run_id, e.returncode)
                with open(err_log, "a") as ef:
                    ef.write(f"{datetime.now().isoformat()}\t{run_id}\t{e.returncode}\n")


def collapse_all_regions(cfg, outdir, dry=False):
    """
    For each region in cfg['regions'], for each sample/corr_run under
      <outdir>/regions/<region_id>/<sample>/<corr_run>/…,
    run exactly the option‐dicts in cfg['collapse_options'], placing outputs under:
      <outdir>/results/<region_id>/<sample>/<corr_run>/<run_id>/…
    """
    genome    = cfg["genome"]
    # default to 4 if threads not provided or zero
    threads   = cfg.get("threads") or 4
    junctions = cfg.get("junctions")
    flag_sets = cfg["collapse_options"]  # list of dicts

    for chrom, spans in cfg["regions"].items():
        for span in spans:
            start, end   = map(int, span.split(":", 1))
            region_id    = f"{chrom}_{start}_{end}"
            slice_root   = os.path.join(outdir, "regions", region_id)

            # walk each sample/corr_run in the slice tree
            pattern = os.path.join(slice_root, "*", "*")
            for sample_run in sorted(glob.glob(pattern)):
                sample   = os.path.basename(os.path.dirname(sample_run))
                corr_run = os.path.basename(sample_run)

                # per‐run config
                rcfg = {
                    "genome":       genome,
                    "threads":      threads,
                    "query_reg":    os.path.join(
                        sample_run, f"{chrom}_{start}-{end}.bed"
                    ),
                    "reads_fa_reg": os.path.join(
                        sample_run, f"{chrom}_{start}-{end}.fasta"
                    ),
                    "gtf_reg":      os.path.join(
                        sample_run, f"{chrom}_{start}-{end}.gtf"
                    ),
                    "junctions_reg": (
                        os.path.join(
                            sample_run, f"{chrom}_{start}-{end}.tab"
                        ) if junctions else None
                    )
                }

                base_cmd     = flair_cmd(rcfg)
                results_base = os.path.join(
                    outdir, "results", region_id, sample, corr_run
                )
                logs_base    = os.path.join(
                    outdir, "logs",    region_id, sample, corr_run
                )

                logging.info(
                    "=== Collapsing %s → %s / %s (threads=%d) ===",
                    region_id, sample, corr_run, threads
                )
                run_jobs(flag_sets, base_cmd, results_base, logs_base, dry)
