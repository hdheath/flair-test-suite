#!/usr/bin/env python3
"""
src/flair_test_suite/align.py

Run FLAIR `align` using shared build/skip/write helpers,
skip runs already completed, and write a completion marker with metadata.
Activates the specified Conda environment in a bash shell for each align.
"""
import logging
import subprocess
from pathlib import Path

from flair_test_suite.build_completion_file import (
    build_metadata,
    should_skip,
    write_marker,
)

def run_align(cfg):
    """
    Perform the alignment stage using FLAIR based on `cfg`.
    Activates the Conda environment `cfg.run.conda_env`, then runs `flair align`.

    Args:
        cfg: Namespace from load_config(), with fields:
            sample_name, input_root, data_dir, run_id,
            run.conda_env, run.steps.align.flags, run.version
    """
    base_out = Path("outputs") / cfg.sample_name / "align"

    input_root = Path(cfg.input_root)
    data = cfg.data_dir
    genome = (input_root / data.genome_fasta).resolve()
    reads_raw = data.reads_fasta
    read_files = (
        [input_root / reads_raw] if isinstance(reads_raw, str)
        else [input_root / rf for rf in reads_raw]
    )
    read_files = [rf.resolve() for rf in read_files]

    flags = vars(cfg.run.steps.align.flags)
    run_id = str(cfg.run_id)
    conda_env = cfg.run.conda_env

    for read_file in read_files:
        sample = cfg.sample_name

        # build metadata and check for skips
        metadata = build_metadata(
            run_id=run_id,
            sample=sample,
            read_file=str(read_file),
            version=cfg.run.version,
            conda_env=conda_env,
            flags=flags,
        )
        match = should_skip(cfg.sample_name, 'align', sample, metadata)
        if match:
            logging.info(
                "[SKIP] align (same settings already ran as run_id=%s) for sample %s",
                match, sample
            )
            continue

        # run FLAIR align
        run_dir = base_out / run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        flair_cmd = ["flair", "align", "-g", str(genome), "-r", str(read_file)]
        for k, v in sorted(flags.items()):
            opt = f"--{k.replace('_', '-') }"
            if isinstance(v, bool):
                if v:
                    flair_cmd.append(opt)
            else:
                flair_cmd.extend([opt, str(v)])
        flair_cmd.extend(["-o", run_id])
        joined = " ".join(flair_cmd)
        cmd = ["bash", "-lc", f"conda activate {conda_env} && {joined}"]

        logging.info("[RUN] cd %s && %s", run_dir, ' '.join(cmd))
        try:
            subprocess.run(cmd, cwd=run_dir, check=True)
        except FileNotFoundError as e:
            logging.error("Executable not found: %s", e)
            raise RuntimeError(f"Conda or FLAIR not found for env '{conda_env}'") from e
        except subprocess.CalledProcessError as e:
            logging.error(
                "FLAIR align failed for run_id=%s sample=%s: %s",
                run_id, sample, e
            )
            continue

        # write marker
        write_marker(cfg.sample_name, 'align', metadata)
