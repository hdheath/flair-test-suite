import types
import sys
from pathlib import Path

import pytest

sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

sys.modules.setdefault("intervaltree", types.SimpleNamespace(IntervalTree=object))
sys.modules.setdefault("pysam", types.SimpleNamespace(AlignedSegment=object))

from flair_test_suite.config_schema import Config as Cfg
from flair_test_suite.qc import ted


def _make_cfg(tmp_path: Path, with_peaks: bool = True) -> Cfg:
    data_dir = tmp_path / "test_data"
    data_dir.mkdir()
    for name in ["exp5.bed", "exp3.bed", "ref5.bed", "ref3.bed"]:
        (data_dir / name).write_text("chr1\t0\t1\t.\t0\t+\n")

    run_cfg = {
        "version": "1",
        "conda_env": "env",
        "work_dir": "outputs",
        "data_dir": "test_data",
        "stages": [{"name": "collapse", "requires": []}],
    }
    if with_peaks:
        run_cfg.update(
            {
                "experiment_5_prime_regions_bed_file": "exp5.bed",
                "experiment_3_prime_regions_bed_file": "exp3.bed",
                "reference_5_prime_regions_bed_file": "ref5.bed",
                "reference_3_prime_regions_bed_file": "ref3.bed",
            }
        )

    cfg_dict = {"run_id": "run1", "run": run_cfg, "qc": {"collapse": {"TED": {"window": 10}}}}
    return Cfg.model_validate(cfg_dict)


def test_collect_uses_run_level_peaks(tmp_path: Path, monkeypatch):
    cfg = _make_cfg(tmp_path, with_peaks=True)
    repo_root = tmp_path
    run_id = "run1"
    stage_uuid = "123abc"
    stage_dir = repo_root / "outputs" / run_id / "collapse" / stage_uuid
    stage_dir.mkdir(parents=True)
    iso_bed = stage_dir / f"{run_id}.isoforms.bed"
    iso_bed.write_text("chr1\t100\t200\tread_ENSG000001.1\t0\t+\n")

    captured = {}

    def fake_metrics(iso_bed_path, peaks, window, audit_rows, tag_ctx):
        captured.update(peaks)
        return {}

    monkeypatch.setattr(ted, "_tss_tts_metrics_full", fake_metrics)

    out_dir = stage_dir / "qc" / "ted"
    ted.collect(stage_dir, cfg, out_dir=out_dir)

    data_dir = repo_root / "test_data"
    assert captured["prime5"] == data_dir / "exp5.bed"
    assert captured["prime3"] == data_dir / "exp3.bed"
    assert captured["ref_prime5"] == data_dir / "ref5.bed"
    assert captured["ref_prime3"] == data_dir / "ref3.bed"


def test_collect_missing_run_level_peaks_errors(tmp_path: Path, monkeypatch):
    cfg = _make_cfg(tmp_path, with_peaks=False)
    repo_root = tmp_path
    run_id = "run1"
    stage_uuid = "123abc"
    stage_dir = repo_root / "outputs" / run_id / "collapse" / stage_uuid
    stage_dir.mkdir(parents=True)
    iso_bed = stage_dir / f"{run_id}.isoforms.bed"
    iso_bed.write_text("chr1\t100\t200\tread_ENSG000001.1\t0\t+\n")

    with pytest.raises(RuntimeError):
        ted.collect(stage_dir, cfg, out_dir=stage_dir / "qc" / "ted")
