# src/flair_test_suite/core/paths.py
# -----------------------------------
# Utilities for constructing and managing file paths for each pipeline stage.
# PathBuilder encapsulates the logic for where stage outputs, markers, and
# sidecar files should live on disk.

from __future__ import annotations
import hashlib  # for checksum helper
from pathlib import Path  # for filesystem path operations
import os 


class PathBuilder:
    """
    Builds canonical directories and filenames for each stage run, based on:
      - work_dir: root directory for all outputs
      - sample:   sample name string
      - stage:    stage name (e.g., 'align', 'correct')
      - signature: unique hash identifying this particular run

    Also provides helper methods for writing markers and computing checksums.

    The `metadata` attribute (populated at runtime by StageBase) can store
    arbitrary key/value pairs, such as QC metrics, for downstream stages.
    """

    # Placeholder for per-run metadata (populated dynamically):
    metadata: dict[str, any] = {}

    def __init__(
        self,
        work_dir: Path,
        sample: str,
        stage: str,
        signature: str,
    ):
        # Base output directory, expanded (~) and resolved to an absolute path
        self.base = Path(work_dir).expanduser().resolve()
        # Sample identifier (used as a sub-directory)
        self.sample = sample
        # Stage identifier (used as a sub-directory)
        self.stage = stage
        # Unique signature (used as the final sub-directory)
        self.signature = signature

    @property
    def stage_dir(self) -> Path:
        """
        Directory for all outputs of this stage/run:
        <base>/<sample>/<stage>/<signature>
        """
        return self.base / self.sample / self.stage / self.signature

    def out(self, *parts: str | Path) -> Path:
        """
        Helper to construct a path under the stage_dir.
        Usage: pb.out("subdir", "file.txt") -> <stage_dir>/subdir/file.txt
        """
        return self.stage_dir.joinpath(*map(str, parts))

    def marker(self) -> Path:
        """
        Path to the JSON completion marker file for this stage.
        """
        return self.stage_dir / ".completed.json"

    @staticmethod
    def sha256(fp: Path, bs: int = 2**20) -> str:
        """
        Compute the SHA-256 checksum of a file, reading in binary chunks.

        Parameters:
          - fp: Path to the file to hash
          - bs: block size (bytes) to read at a time (default 1 MiB)

        Returns:
          - hex digest string
        """
        h = hashlib.sha256()
        with open(fp, "rb") as f:
            # read in blocks until EOF
            for chunk in iter(lambda: f.read(bs), b""):
                h.update(chunk)
        return h.hexdigest()

    @staticmethod
    def sha256_str(txt: str) -> str:
        """
        Compute the SHA-256 checksum of a UTF-8 string.

        Parameters:
          - txt: input string

        Returns:
          - hex digest string
        """
        return hashlib.sha256(txt.encode()).hexdigest()
