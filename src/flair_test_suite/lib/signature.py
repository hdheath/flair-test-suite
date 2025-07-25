# src/flair_test_suite/core/signature.py
# ---------------------------------------
# Utilities to compute a concise signature for a stage run,
# write and load JSON completion markers, and check completeness.

from __future__ import annotations

import json           # for reading/writing JSON metadata
import hashlib        # for computing SHA-1 signatures
import base64         # optionally available for encoding, not used here
from datetime import datetime, timezone  # for timestamping markers
import platform       # to record host information in markers
from pathlib import Path  # for filesystem path operations

# Import PathBuilder to know where to write marker files
from .paths import PathBuilder


def compute_signature(
    tool_version: str,
    flags: str,
    hashes: list[str]
) -> str:
    """
    Generate a stable, short signature string for a stage execution.

    Combines:
      - tool_version: version string of the external tool (e.g., "flair 3.0.0")
      - flags:       concatenated CLI flags used
      - hashes:      list of input file or scalar hashes

    Steps:
      1. Deduplicate and sort `hashes` for stable ordering.
      2. Create a payload string: "tool_version|flags|hash1|hash2|..."
      3. Compute SHA-1 digest of the payload and take first 12 hex chars.

    Returns:
      - 12-character hex string unique to this combination.
    """
    # Remove duplicates and sort for deterministic signature
    uniq = sorted(set(hashes))
    # Concatenate all parts into a single string
    payload = f"{tool_version}|{flags}|{'|'.join(uniq)}"
    # Compute SHA-1 and truncate to 12 characters
    return hashlib.sha1(payload.encode()).hexdigest()[:12]


def write_marker(pb: PathBuilder, meta: dict):
    """
    Write the completion metadata `meta` to a JSON file in the stage directory.

    - Ensures the directory exists.
    - Writes to `.completed.json` with 2-space indentation.

    Parameters:
      - pb:   PathBuilder instance identifying the stage run.
      - meta: dict of metadata (including QC, runtime, etc.).
    """
    # Ensure output directory exists
    pb.stage_dir.mkdir(parents=True, exist_ok=True)
    # Write metadata JSON to the marker file
    with open(pb.marker(), "w") as fh:
        json.dump(meta, fh, indent=2)


def qc_sidecar_path(stage_dir: Path, stage_name: str) -> Path:
    """
    Construct the path to the QC sidecar TSV for a given stage.

    Parameters:
      - stage_dir: base directory for the stage signature
      - stage_name: name of the stage (used as filename prefix)

    Returns:
      - <stage_dir>/<stage_name>_qc.tsv
    """
    return stage_dir / f"{stage_name}_qc.tsv"


def load_marker(stage_dir: Path) -> dict | None:
    """
    Read the JSON completion marker if it exists.

    Parameters:
      - stage_dir: directory containing `.completed.json`

    Returns:
      - Parsed dict if marker exists, else None.
    """
    m = stage_dir / ".completed.json"
    if m.exists():
        return json.loads(m.read_text())
    return None


def is_complete(
    stage_dir: Path,
    outputs: list[Path],
    needs_qc: bool
) -> bool:
    """
    Check whether a stage is fully complete:
      - `.completed.json` exists
      - All expected output files exist
      - If QC is needed: QC TSV exists and marker contains non-empty 'qc' block

    Parameters:
      - stage_dir: where the stage outputs live
      - outputs:   list of Path objects expected
      - needs_qc:  whether this stage has a QC collector

    Returns:
      - True if all conditions are met, False otherwise.
    """
    # 1) Marker must exist
    marker = stage_dir / ".completed.json"
    if not marker.exists():
        return False
    # 2) Each expected output file must exist
    if not all(p.exists() for p in outputs):
        return False
    # 3) If QC is required, ensure TSV and QC block are present
    if needs_qc:
        qc_tsv = stage_dir / f"{stage_dir.name}_qc.tsv"
        if not qc_tsv.exists():
            return False
        # Load marker JSON and verify 'qc' key has content
        if not json.loads(marker.read_text()).get("qc"):
            return False
    return True
