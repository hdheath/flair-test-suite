#!/usr/bin/env python3
"""
src/flair_test_suite/build_completion_file.py

Helper utilities to assemble metadata, determine if a stage should be skipped,
and write a completion marker for any pipeline stage.
"""
import logging
from pathlib import Path

from flair_test_suite.check_completion_markers import find_matching_run, marker_path


def build_metadata(run_id: str,
                   sample: str,
                   read_file: str,
                   version: str,
                   conda_env: str,
                   flags: dict) -> dict:
    """
    Assemble the metadata dict for any stage.

    Args:
        run_id: the run identifier (stringified)
        sample: the sample name
        read_file: absolute path to the input read file
        version: the pipeline/software version
        conda_env: the conda environment used
        flags: dict of flag_name->flag_value

    Returns:
        A dict mapping metadata keys to their string values.
    """
    # canonicalize flags: drop False/None/"", sort, join
    parts = [f"{k}={v}" for k, v in sorted(flags.items()) if v not in (False, None, "")]
    flags_str = "; ".join(parts)

    return {
        "Config Run ID": run_id,
        "Sample": sample,
        "Read file": read_file,
        "Version": version,
        "Conda Env": conda_env,
        "Flags": flags_str,
    }


def should_skip(sample_name: str,
                stage: str,
                sample: str,
                metadata: dict) -> str | None:
    """
    Determine whether an identical run (ignoring run ID) already completed.

    Returns the matching prior run_id if found, else None.
    """
    return find_matching_run(sample_name, stage, sample, metadata)


def write_marker(sample_name: str,
                 stage: str,
                 metadata: dict) -> Path:
    """
    Write a completion marker file based on metadata.

    Returns the Path to the written marker.
    """
    run_id = metadata["Config Run ID"]
    sample = metadata.get("Sample")
    marker = marker_path(sample_name, stage, sample, run_id)
    marker.parent.mkdir(parents=True, exist_ok=True)

    # one key:value per line
    lines = [f"{k}: {v}" for k, v in metadata.items()]
    marker.write_text("\n".join(lines))
    logging.info("[DONE] %s complete; marker at %s", stage, marker)
    return marker
