# src/flair_test_suite/core/input_hash.py
# ------------------------------------------------
# Utilities for hashing files and scalar values to
# produce stable input signatures for pipeline stages.

from __future__ import annotations
import hashlib  # for SHA-256 hashing
import json     # to serialize scalars consistently
from pathlib import Path  # for filesystem paths
from typing import Iterable, List, Union  # type annotations

# Explicitly export the public functions of this module
__all__ = ["hash_path", "hash_scalar", "hash_many"]


def _sha256_bytes(data: bytes) -> str:
    """
    Compute the SHA-256 hex digest of raw bytes.

    Parameters:
      - data: raw bytes to hash

    Returns:
      - lowercase hex string of the SHA-256 digest
    """
    return hashlib.sha256(data).hexdigest()


def hash_path(p: Path) -> str:
    """
    Hash the entire contents of a file using SHA-256.

    Parameters:
      - p: Path to the file (must exist)

    Returns:
      - hex digest string of the file's contents

    Raises:
      - any IO error if the file is missing or unreadable
    """
    # Read raw bytes from the file and hash them
    return _sha256_bytes(p.read_bytes())


def hash_scalar(x: Union[str, int, float, bool]) -> str:
    """
    Hash a scalar value by JSON-encoding it, ensuring a stable
    representation across Python runs.

    Parameters:
      - x: a primitive type (str, int, float, or bool)

    Returns:
      - hex digest string of the JSON-encoded bytes
    """
    # json.dumps with sort_keys keeps output consistent for dict-like inputs,
    # though here x is a scalar; encode to UTF-8 before hashing
    return _sha256_bytes(json.dumps(x, sort_keys=True).encode())


def hash_many(
    items: Iterable[Union[Path, str, int, float, bool]]
) -> List[str]:
    """
    Compute hashes for a mixed iterable of file Paths and scalar values.
    Returns a sorted, de-duplicated list of hex digest strings.

    Parameters:
      - items: iterable containing Path objects or primitive scalars

    Returns:
      - sorted list of unique hash strings
    """
    out: set[str] = set()
    for it in items:
        if isinstance(it, Path):
            # For Path, resolve to absolute path then hash file contents
            out.add(hash_path(it.resolve()))
        else:
            # For scalars (or raw strings), hash their JSON representation
            out.add(hash_scalar(it))
    # sorted() ensures a deterministic order for signature computation
    return sorted(out)

