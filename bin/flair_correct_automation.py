#!/usr/bin/env python3
"""
flair_automate/flair_correct_automation.py

Run FLAIR `correct` for each `flair align` output and each set of
user-specified flags, but skip any runs whose primary output already exists.
"""
import os
import subprocess
import logging

def correct_all(cfg, outdir, dry_run=False):
    """
    For each sample in cfg['reads_fasta'], for each align‐run under
    <outdir>/align/<sample>/<align_run>/, and for each dict in
    cfg.get('correct_options', []), run:

        flair correct \
          -q <align_dir>/<sample>.bed \
          -f <cfg['gtf']> \
          -j <cfg['junctions']> \
          -g <cfg['genome']> \
          [flags...] \
          -o <sample>

    in its own subdirectory under
    <outdir>/correct/<sample>/<align_run>/<run_id>/,
    skipping if `<sample>.corrected.bed` already exists.
    """
    genome       = cfg["genome"]
    gtf          = cfg["gtf"]
    sj           = cfg.get("junctions")
    correct_sets = cfg.get("correct_options", [])

    for read in cfg["reads_fasta"]:
        sample     = os.path.splitext(os.path.basename(read))[0]
        align_base = os.path.join(outdir, "align", sample)
        if not os.path.isdir(align_base):
            logging.warning("No align directory for sample %s, skipping correct", sample)
            continue

        for align_run in sorted(os.listdir(align_base)):
            align_dir = os.path.join(align_base, align_run)
            bed       = os.path.join(align_dir, f"{sample}.bed")
            if not os.path.exists(bed):
                logging.warning("  Missing BED in %s, skipping", align_dir)
                continue

            for opts in correct_sets:
                # build run_id from flags, skipping --threads
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

                corr_dir     = os.path.join(outdir, "correct", sample, align_run, run_id)
                expected_bed = os.path.join(corr_dir, f"{sample}.corrected.bed")

                # 1) Skip if already done
                if not dry_run and os.path.exists(expected_bed):
                    logging.info("  [SKIP] %s/%s (found %s)", sample, run_id, expected_bed)
                    continue

                # 2) Prepare output directory
                os.makedirs(corr_dir, exist_ok=True)

                # 3) Build command (including threads if present)
                cmd = ["flair", "correct",
                       "-q", bed,
                       "-f", gtf]
                if sj:
                    cmd += ["-j", sj]
                cmd += ["-g", genome]
                for flag in sorted(opts):
                    val = opts[flag]
                    if isinstance(val, bool):
                        if val:
                            cmd.append(flag)
                    else:
                        cmd += [flag, str(val)]
                cmd += ["-o", sample]

                # 4) Log & execute
                if dry_run:
                    logging.info("  [DRY] cd %s && %s", corr_dir, " ".join(cmd))
                else:
                    logging.info("  [RUN] cd %s && %s", corr_dir, " ".join(cmd))
                    subprocess.run(cmd, cwd=corr_dir, check=True)
