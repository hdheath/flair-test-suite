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

def _normalize_metadata(meta: dict) -> dict:
    """
    Return a new dict with:
      - "Config Run ID" removed
      - Flags stripped of any threads=... and the remainder sorted & joined
    """
    norm = {}
    for key, val in meta.items():
        if key == "Config Run ID":
            continue
        if key == "Flags":
            # split on ';', strip whitespace, drop threads=, sort the rest
            parts = [p.strip() for p in val.split(";") if p.strip() and not p.strip().startswith("threads=")]
            parts.sort()
            norm["Flags"] = ";".join(parts)
        else:
            norm[key] = val
    return norm

def find_matching_run(sample_name: str,
                      stage: str,
                      sample: str,
                      current_metadata: dict) -> str | None:
    """
    Scan outputs/<sample_name>/<stage>/<run_id> for any marker whose
    normalized metadata (Sample, Read file, Version, Conda Env, Flags)
    exactly equals the normalized current_metadata.  Return that run_id,
    or None if there is no match.
    """
    base = Path("outputs") / sample_name / stage
    if not base.exists():
        return None

    target = _normalize_metadata(current_metadata)

    for candidate in base.iterdir():
        if not candidate.is_dir():
            continue
        run_id = candidate.name

        marker = base / run_id / f"{run_id}.{stage}.completed.txt"
        if not marker.exists():
            # skip silently if no marker
            continue

        old = parse_marker(sample_name, stage, sample, run_id)
        old_norm = _normalize_metadata(old)

        if old_norm == target:
            logging.info("Found matching run %s for %s/%s", run_id, sample_name, stage)
            return run_id

    return None