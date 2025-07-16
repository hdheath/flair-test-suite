#!/usr/bin/env python3
"""
src/flair_test_suite/align.py

Run FLAIR `align` for each input read file using flags from config,
skip runs already completed, and write a completion marker with metadata.
This uses `conda activate <env>` inside a bash login shell so the specified
Conda environment is activated for each align command.
"""
import logging
import subprocess
from pathlib import Path

from flair_test_suite.check_completion_markers import is_run_completed, marker_path


def run_align(cfg, dry_run=False):
    """
    Perform the alignment stage using FLAIR based on `cfg`.
    Activates the Conda environment `cfg.run.conda_env` in a bash login shell,
    then runs the `flair align` command.

    Args:
        cfg: Namespace from load_config(), with fields:
            sample_name, input_root, data_dir, run_id,
            run.conda_env, run.steps.align.flags
        dry_run (bool): If True, prints commands but does not execute.
    """
    base_out = Path("outputs") / cfg.sample_name / "align"
    input_root = Path(cfg.input_root)
    data = cfg.data_dir
    genome = input_root / data.genome_fasta
    reads = data.reads_fasta
    read_files = ([input_root / reads]
                  if isinstance(reads, str)
                  else [input_root / rf for rf in reads])

    flags = vars(cfg.run.steps.align.flags)
    run_id = str(cfg.run_id)
    conda_env = cfg.run.conda_env

    for read_file in read_files:
        sample = read_file.stem
        run_dir = base_out / sample / run_id
        marker = marker_path(cfg.sample_name, 'align', sample, run_id)

        if not dry_run and is_run_completed(cfg.sample_name, 'align', sample, run_id):
            logging.info(f"[SKIP] align {run_id} for sample {sample} (marker exists)")
            continue

        run_dir.mkdir(parents=True, exist_ok=True)

        # Build the bash -lc command to activate env then run flair
        flair_cmd = ["flair", "align", "-g", str(genome), "-r", str(read_file)]
        for k, v in sorted(flags.items()):
            opt = f"--{k.replace('_','-')}"
            if isinstance(v, bool):
                if v:
                    flair_cmd.append(opt)
            else:
                flair_cmd.extend([opt, str(v)])
        flair_cmd.extend(["-o", run_id])
        joined = " ".join(flair_cmd)
        cmd = ["bash", "-lc", f"conda activate {conda_env} && {joined}"]

        if dry_run:
            logging.info(f"[DRY] cd {run_dir} && {' '.join(cmd)}")
            continue

        logging.info(f"[RUN] cd {run_dir} && {' '.join(cmd)}")
        try:
            subprocess.run(cmd, cwd=run_dir, check=True)
        except FileNotFoundError as e:
            logging.error(f"Executable not found: {e}")
            raise RuntimeError(f"Conda or flair not found for env '{conda_env}'") from e
        except subprocess.CalledProcessError as e:
            logging.error(f"FLAIR align failed for run_id={run_id} sample={sample}: {e}")
            continue

        # Write completion marker
        marker.parent.mkdir(parents=True, exist_ok=True)
        try:
            content = [
                f"Config Run ID: {run_id}",
                f"Sample: {sample}",
                f"Read file: {read_file}",
                "Flags: " + "; ".join(f"{k}={v}" for k, v in flags.items())
            ]
            marker.write_text("\n".join(content))
            logging.info(f"[DONE] align complete; marker at {marker}")
        except Exception as e:
            logging.error(f"Failed writing marker {marker}: {e}")
