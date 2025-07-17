#!/usr/bin/env python3
"""
src/flair_test_suite/stages/align.py

Run the FLAIR "align" stage using shared build/skip/write helpers.
This module exposes a standardized `run(flags, cfg, prev_run_id)` entrypoint
so it can be dynamically invoked by the CLI dispatcher and chain run_ids.
"""
import logging
import subprocess
from pathlib import Path

from flair_test_suite.build_completion_file import (
    build_metadata,
    should_skip,
    write_marker,
)


def run(flags, cfg, prev_run_id=None):
    """
    Perform the alignment stage using FLAIR based on `cfg` and `flags`.

    Args:
        flags:        SimpleNamespace of CLI flags from `run.stages.flags`.
        cfg:          Namespace from load_config(), with fields:
                      sample_name, input_root, data_dir, run_id,
                      run.version, run.conda_env, etc.
        prev_run_id:  If provided, the run_id to inspect/skip before using cfg.run_id.

    Returns:
        str: the run_id whose .completed marker was used or created.
    """
    # Determine which run_id to use
    run_id = prev_run_id or str(cfg.run_id)

    # Prepare directories and inputs
    base_out = Path("outputs") / cfg.sample_name / "align"
    input_root = Path(cfg.input_root)
    data = cfg.data_dir
    genome = (input_root / data.genome_fasta).resolve()

    # Determine read files (single or list)
    reads_raw = data.reads_fasta
    if isinstance(reads_raw, str):
        read_files = [input_root / reads_raw]
    else:
        read_files = [input_root / rf for rf in reads_raw]
    read_files = [rf.resolve() for rf in read_files]

    # Convert flags namespace to dict and extract optional aligner override
    flag_dict = vars(flags)
    # Auto-resolve boolean flags to actual file paths
    if flag_dict.get("junction_bed") is True:
        flag_dict["junction_bed"] = str((input_root / data.junctions).resolve())
    aligner = flag_dict.pop('aligner', None)
    cmd_base = aligner.split() if aligner else ["flair", "align"]
    conda_env = cfg.run.conda_env

    for read_file in read_files:
        sample = cfg.sample_name

        # Build metadata and check skip logic
        metadata = build_metadata(
            run_id=run_id,
            sample=sample,
            read_file=str(read_file),
            version=cfg.run.version,
            conda_env=conda_env,
            flags=flag_dict,
        )
        skip_run = should_skip(cfg.sample_name, 'align', sample, metadata)
        if skip_run:
            logging.info(
                "[SKIP] align (using existing run_id=%s) for sample %s",
                skip_run, sample
            )
            return str(skip_run)

        # Create run directory
        run_dir = base_out / run_id
        run_dir.mkdir(parents=True, exist_ok=True)

        # Build the FLAIR align command
        flair_cmd = cmd_base + ["-g", str(genome), "-r", str(read_file)]
        # Forward all flags, including junction_bed
        for k, v in sorted(flag_dict.items()):
            opt = f"--{k}"
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

        # Write completion marker for new run
        write_marker(cfg.sample_name, 'align', metadata)
        return run_id

    # No reads processed? return whatever run_id we have
    return run_id
