#!/usr/bin/env python3
"""
flair_automate/flair_align_automation.py

Run FLAIR `align` for each input read file and each set of user-specified flags,
but skip any runs whose primary output already exists.
"""
import os
import subprocess
import logging

def align_all(cfg, outdir, dry_run=False):
    """
    For each file in cfg['reads_fasta'] and each dict in cfg['align_options'],
    run:
        flair align -g <genome> -r <reads> [flags...] -o <prefix>
    in its own subdirectory under <outdir>/align/<sample>/<run_id>/,
    skipping if <run_dir>/<sample>.bam already exists.
    """
    genome     = cfg["genome"]
    reads      = cfg["reads_fasta"]
    align_sets = cfg.get("align_options", [])

    for read in reads:
        sample = os.path.splitext(os.path.basename(read))[0]

        for opts in align_sets:
            # build run_id from flags, but skip --threads
            parts = []
            for flag in sorted(opts):
                if flag == "--threads":
                    continue
                val = opts[flag]
                if isinstance(val, bool):
                    if val:
                        parts.append(flag.lstrip("-"))
                else:
                    parts.append(f"{flag.lstrip('-')}-{val}")
            run_id = "_".join(parts) if parts else "default"

            run_dir      = os.path.join(outdir, "align", sample, run_id)
            expected_bam = os.path.join(run_dir, f"{sample}.bam")

            # SKIP if we've already aligned
            if not dry_run and os.path.exists(expected_bam):
                logging.info("  [SKIP] %s (found %s)", run_id, expected_bam)
                continue

            os.makedirs(run_dir, exist_ok=True)

            # build command (including threads)
            cmd = ["flair", "align", "-g", genome, "-r", read]
            for flag in sorted(opts):
                val = opts[flag]
                if isinstance(val, bool):
                    if val:
                        cmd.append(flag)
                else:
                    cmd.extend([flag, str(val)])
            cmd += ["-o", sample]

            # log & run
            if dry_run:
                logging.info("  [DRY] cd %s && %s", run_dir, " ".join(cmd))
            else:
                logging.info("  [RUN] cd %s && %s", run_dir, " ".join(cmd))
                subprocess.run(cmd, cwd=run_dir, check=True)
