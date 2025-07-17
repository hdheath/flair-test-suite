"""
Tier‑1 smoke test for the Align stage.

Checks:
  • .completed.json exists
  • exit_code == 0
  • runtime < cfg.qc.align.max_runtime_hours
  • mapped_pct >= cfg.qc.align.min_mapped_pct
  • retained_pct >= cfg.qc.align.min_retained_pct
  • QC side‑car files (bam, bed, tsv) are present
"""
from __future__ import annotations
import json
from datetime import timedelta
from pathlib import Path

import pytest
from flair_test_suite.config_loader import load_config

CFG_PATH = Path(__file__).parents[1] / "config" / "human_wtc11_align.toml"
SAMPLES  = ["WTC11_demo"]          # extend with more sample names as you add data


# ---------- fixtures -------------------------------------------------
@pytest.fixture(scope="session")
def cfg():
    return load_config(CFG_PATH)


def _stage_dir(cfg, sample: str):
    root = Path(cfg.run.work_dir)
    align_root = root / sample / "align"
    # assume exactly one signature folder for this sample
    sig_dirs = list(align_root.iterdir())
    assert sig_dirs, f"No align signature folder for {sample}"
    return sig_dirs[0]


def _load_meta(stage_dir: Path):
    marker = stage_dir / ".completed.json"
    assert marker.exists(), f"Missing .completed.json in {stage_dir}"
    return json.loads(marker.read_text())


# ---------- parametrised test ---------------------------------------
@pytest.mark.parametrize("sample", SAMPLES)
def test_align_qc(cfg, sample):
    stage_dir = _stage_dir(cfg, sample)
    meta      = _load_meta(stage_dir)

    # basic sanity
    assert meta["exit_code"] == 0, "Align exited with non‑zero status"

    max_rt = timedelta(hours=cfg.qc.align.max_runtime_hours).total_seconds()
    assert meta["runtime_sec"] < max_rt, (
        f"Align runtime {meta['runtime_sec']}s exceeds {max_rt}s"
    )

    qc = meta.get("qc", {})
    assert qc, "QC metrics missing from marker!"

    mpct  = qc.get("mapped_pct")
    r_pct = qc.get("retained_pct")
    assert mpct is not None, "mapped_pct missing"
    assert r_pct is not None, "retained_pct missing"

    assert mpct >= cfg.qc.align.min_mapped_pct, (
        f"mapped_pct {mpct}% < threshold {cfg.qc.align.min_mapped_pct}%"
    )
    assert r_pct >= cfg.qc.align.min_retained_pct, (
        f"retained_pct {r_pct}% < threshold {cfg.qc.align.min_retained_pct}%"
    )

    # outputs present
    bam = stage_dir / f"{sample}_flair.bam"
    bed = stage_dir / f"{sample}_flair.bed"
    tsv = stage_dir / "align_qc.tsv"
    for fp in (bam, bed, tsv):
        assert fp.exists(), f"Expected output missing: {fp}"

