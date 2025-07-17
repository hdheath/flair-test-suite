# tests/align/test_align_qc.py
from datetime import timedelta
from pathlib import Path

def test_align_qc(cfg, sample, stage_dir, align_meta):
    meta = align_meta
    assert meta["exit_code"] == 0

    max_rt = timedelta(hours=cfg.qc.align.max_runtime_hours).total_seconds()
    assert meta["runtime_sec"] < max_rt

    qc = meta["qc"]
    assert qc["mapped_pct"] >= cfg.qc.align.min_mapped_pct
    assert qc["retained_pct"] >= cfg.qc.align.min_retained_pct

    # outputs present
    for fname in [f"{sample}_flair.bam",
                  f"{sample}_flair.bed",
                  "align_qc.tsv"]:
        assert (stage_dir / fname).exists(), f"Missing {fname}"

