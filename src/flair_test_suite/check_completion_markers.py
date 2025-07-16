"""
src/flair_test_suite/check_completion_markers.py

Common utilities for checking and parsing run completion markers.
"""
from pathlib import Path
import logging

def marker_path(sample_name: str, stage: str, sample: str, run_id: str) -> Path:
    """
    Construct the expected marker file path for a given stage/run.
    e.g. outputs/<sample_name>/<stage>/<sample>/<run_id>/<run_id>.<stage>.completed.txt
    """
    return Path("outputs") / sample_name / stage / sample / run_id / f"{run_id}.{stage}.completed.txt"


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
