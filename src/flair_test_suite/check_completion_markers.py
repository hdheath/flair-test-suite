"""
src/flair_test_suite/check_completion_markers.py

Common utilities for checking and parsing run completion markers.
"""
from pathlib import Path
import logging

def marker_path(sample_name: str, stage: str, sample: str, run_id: str) -> Path:
    """
    Construct the expected marker file path for a given stage/run.
    e.g. outputs/<sample_name>/<stage>/<run_id>/<run_id>.<stage>.completed.txt
    """
    return (Path("outputs") 
        / sample_name 
        / stage 
        / run_id 
        / f"{run_id}.{stage}.completed.txt")


def is_run_completed(sample_name: str, stage: str, sample: str, run_id: str) -> bool:
    """
    Check if the completion marker for this stage/run exists.
    """
    marker = marker_path(sample_name, stage, sample, run_id)
    exists = marker.exists()
    logging.debug("Checking marker %s exists: %s", marker, exists)
    return exists


def parse_marker(sample_name: str, stage: str, sample: str, run_id: str) -> dict:
    """
    Read and parse the completion marker file into a dict of metadata.
    Returns an empty dict if marker missing.
    """
    marker = marker_path(sample_name, stage, sample, run_id)
    data = {}
    if not marker.exists():
        logging.warning("Marker file not found: %s", marker)
        return data
    for line in marker.read_text().splitlines():
        if ':' in line:
            key, val = line.split(':', 1)
            data[key.strip()] = val.strip()
    return data


def _sanitize_flags(flags_str: str) -> str:
    """
    Remove any 'threads=...' entries from a semicolon‑separated Flags string.
    """
    parts = [p for p in flags_str.split(';') if not p.strip().startswith("threads=")]
    return ";".join(p.strip() for p in parts if p.strip())

def find_matching_run(sample_name: str,
                      stage: str,
                      sample: str,
                      current_metadata: dict) -> str | None:
    """
    Scan all run_id subdirs under outputs/<sample_name>/<stage>/,
    parse each marker’s metadata, and if one matches current_metadata
    (ignoring run-id and any threads=... in Flags), return that prior run_id.
    Otherwise return None.
    """
    base = Path("outputs") / sample_name / stage
    if not base.exists():
        return None

    # Prepare the target metadata, minus run-id, and with sanitized Flags
    target = {}
    for k, v in current_metadata.items():
        if k == "Config Run ID":
            continue
        if k == "Flags":
            target["Flags"] = _sanitize_flags(v)
        else:
            target[k] = v

    for candidate in base.iterdir():
        if not candidate.is_dir():
            continue
        run_id = candidate.name
        old = parse_marker(sample_name, stage, sample, run_id)
        if not old:
            continue

        # Sanitize the old metadata the same way
        old_sanitized = {}
        for k, v in old.items():
            if k == "Config Run ID":
                continue
            if k == "Flags":
                old_sanitized["Flags"] = _sanitize_flags(v)
            else:
                old_sanitized[k] = v

        # Compare sanitized old → target
        if old_sanitized == target:
            logging.info("Found matching run %s (ignoring threads)", run_id)
            return run_id

    return None