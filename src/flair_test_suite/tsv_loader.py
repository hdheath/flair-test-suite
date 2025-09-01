from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, List

from .config_schema import Config


def _coerce_scalar(val: str) -> Any:
    s = val.strip()
    if s.lower() in ("true", "false"):
        return s.lower() == "true"
    try:
        if s.startswith("0") and s != "0":
            raise ValueError
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        return val


def _split_flags_cell(s: str) -> List[str]:
    s = (s or "").strip()
    if not s:
        return []
    reader = csv.reader([s], delimiter=",", quotechar='"', skipinitialspace=True)
    row = next(reader)
    return [tok.strip() for tok in row if tok and tok.strip()]


def build_configs_from_two_section(tsv_path: Path) -> list[Config]:
    """Parse a two-section TSV into Config objects.

    Format:
      key<TAB>value           (global inputs for the test set)
      ...
      [blank or comments]
      stage<TAB>flags         (one row per stage; a row with stage=align starts a new case)
      ...
    """
    lines = tsv_path.read_text().splitlines()
    # Locate headers
    key_header = None
    stage_header = None
    for i, raw in enumerate(lines):
        ln = raw.strip().lower()
        if ln == "key\tvalue":
            key_header = i
        if ln == "stage\tflags":
            stage_header = i
    if key_header is None or stage_header is None or stage_header <= key_header:
        raise ValueError("Two-section TSV requires 'key\tvalue' then 'stage\tflags' headers")

    # Parse inputs
    inputs: Dict[str, Any] = {}
    for raw in lines[key_header + 1: stage_header]:
        if not raw.strip() or raw.lstrip().startswith('#'):
            continue
        parts = raw.split('\t')
        if len(parts) < 2:
            continue
        k = parts[0].strip()
        v = parts[1].strip()
        inputs[k] = _coerce_scalar(v)
    # Basic required keys (environment validated later via model validator)
    required = ["version", "data_dir", "reads_file"]
    missing = [k for k in required if not str(inputs.get(k, "")).strip()]
    if missing:
        raise ValueError("Missing required inputs: " + ", ".join(missing))

    test_set_id = str(inputs.get("test_set_id", "")).strip() or None

    # Parse stages
    def finalize_case(stages_accum: list[dict], out: list[Config]):
        if not stages_accum:
            return
        data = {
            "test_set_id": test_set_id,
            "run": {k: v for k, v in inputs.items() if k != "test_set_id"} | {"stages": stages_accum},
            "qc": {},
        }
        out.append(Config.parse_obj(data))

    cfgs: list[Config] = []
    stages_accum: list[dict] = []
    for raw in lines[stage_header + 1:]:
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        parts = raw.split('\t')
        stage = (parts[0] if parts else "").strip().lower()
        if not stage:
            continue
        if stage == "regionalize":
            raise ValueError("Stage 'regionalize' is not supported; use 'region_test'.")
        flags_cell = parts[1].strip() if len(parts) > 1 else ""
        flags_tokens = _split_flags_cell(flags_cell)

        # A new align row starts a new case
        if stage == "align" and stages_accum:
            finalize_case(stages_accum, cfgs)
            stages_accum = []
        stages_accum.append({"name": stage, "flags": flags_tokens})

    finalize_case(stages_accum, cfgs)
    if not cfgs:
        raise ValueError("No cases parsed; ensure there is at least one 'align' row under 'stage\tflags'.")
    return cfgs
