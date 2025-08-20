import sys
from pathlib import Path

import pandas as pd
import pytest

# Make package importable
sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

from flair_test_suite.qc import ted


def test_collect_single_resolves_peaks(tmp_path: Path, monkeypatch):
    repo_root = tmp_path
    data_dir = repo_root / "test_data"
    data_dir.mkdir()
    for name in ["exp5.bed", "exp3.bed", "ref5.bed", "ref3.bed"]:
        (data_dir / name).write_text("chr1\t0\t1\t.\t0\t+\n")

    run_id = "run1"
    stage_uuid = "123abc"
    stage_dir = repo_root / "outputs" / run_id / "collapse" / stage_uuid
    stage_dir.mkdir(parents=True)

    iso_bed = stage_dir / f"{run_id}.isoforms.bed"
    iso_bed.write_text("chr1\t100\t200\tread_ENSG000001.1\t0\t+\n")
    map_txt = stage_dir / f"{run_id}.isoform.read.map.txt"
    map_txt.write_text(f"{iso_bed.stem}\tread1,read2\n")

    corr_dir = repo_root / "outputs" / run_id / "correct" / "abcd"
    corr_dir.mkdir(parents=True)
    (corr_dir / f"{run_id}_all_corrected.bed").write_text(
        "chr1\t100\t200\tread_ENSG000001.1\t0\t+\n"
    )

    cfg = {
        "run": {"data_dir": "test_data"},
        "qc": {
            "collapse": {
                "TED": {
                    "experiment_5_prime_regions_bed_file": "exp5.bed",
                    "experiment_3_prime_regions_bed_file": "exp3.bed",
                    "reference_5_prime_regions_bed_file": "ref5.bed",
                    "reference_3_prime_regions_bed_file": "ref3.bed",
                    "window": 10,
                }
            }
        },
    }

    captured = {}

    def fake_metrics(iso_bed_path, peaks, window, audit_rows, tag_ctx):
        captured.update(peaks)
        return {
            "5prime_precision": 0.1,
            "5prime_recall": 0.1,
            "5prime_f1": 0.1,
            "3prime_precision": 0.2,
            "3prime_recall": 0.2,
            "3prime_f1": 0.2,
            "ref5prime_precision": 0.3,
            "ref5prime_recall": 0.3,
            "ref5prime_f1": 0.3,
            "ref3prime_precision": 0.4,
            "ref3prime_recall": 0.4,
            "ref3prime_f1": 0.4,
        }

    monkeypatch.setattr(ted, "_tss_tts_metrics_full", fake_metrics)

    ted.collect(stage_dir, cfg)

    assert captured["prime5"] == data_dir / "exp5.bed"
    assert captured["prime3"] == data_dir / "exp3.bed"
    assert captured["ref_prime5"] == data_dir / "ref5.bed"
    assert captured["ref_prime3"] == data_dir / "ref3.bed"

    df = pd.read_csv(stage_dir / "TED.tsv", sep="\t")
    assert df.loc[0, "5prime_precision"] == 0.1
    assert df.loc[0, "3prime_precision"] == 0.2
    assert df.loc[0, "ref5prime_precision"] == 0.3
    assert df.loc[0, "ref3prime_precision"] == 0.4


def test_collect_single_uses_stage_flags(tmp_path: Path, monkeypatch):
    """If the QC block omits peak paths, fall back to run.stage flags."""
    repo_root = tmp_path
    data_dir = repo_root / "test_data"
    data_dir.mkdir()
    for name in ["exp5.bed", "exp3.bed", "ref5.bed", "ref3.bed"]:
        (data_dir / name).write_text("chr1\t0\t1\t.\t0\t+\n")

    run_id = "run1"
    stage_uuid = "123abc"
    stage_dir = repo_root / "outputs" / run_id / "collapse" / stage_uuid
    stage_dir.mkdir(parents=True)

    iso_bed = stage_dir / f"{run_id}.isoforms.bed"
    iso_bed.write_text("chr1\t100\t200\tread_ENSG000001.1\t0\t+\n")

    cfg = {
        "run": {
            "data_dir": "test_data",
            "stages": {
                "collapse": {
                    "flags": {
                        "experiment_5_prime_regions_bed_file": "exp5.bed",
                        "experiment_3_prime_regions_bed_file": "exp3.bed",
                        "reference_5_prime_regions_bed_file": "ref5.bed",
                        "reference_3_prime_regions_bed_file": "ref3.bed",
                    }
                }
            },
        },
        "qc": {"collapse": {"TED": {"window": 10}}},
    }

    captured = {}

    def fake_metrics(iso_bed_path, peaks, window, audit_rows, tag_ctx):
        captured.update(peaks)
        return {}

    monkeypatch.setattr(ted, "_tss_tts_metrics_full", fake_metrics)

    ted.collect(stage_dir, cfg)

    assert captured["prime5"] == data_dir / "exp5.bed"
    assert captured["prime3"] == data_dir / "exp3.bed"
    assert captured["ref_prime5"] == data_dir / "ref5.bed"
    assert captured["ref_prime3"] == data_dir / "ref3.bed"
