import pytest
from pathlib import Path
import sys

# Make package importable and skip if heavy QC deps missing
sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))
pytest.importorskip("pysam")
pytest.importorskip("matplotlib")

from flair_test_suite import qc


def test_qc_registry_has_expected_collectors():
    expected = {"align", "correct", "regionalize", "ted", "transcriptome"}
    assert expected <= set(qc.QC_REGISTRY)
    for func in qc.QC_REGISTRY.values():
        assert callable(func)


def test_write_metrics_creates_tsv(tmp_path: Path):
    metrics = {"one": 1, "two": "dos"}
    qc.write_metrics(tmp_path, "dummy", metrics)
    out = tmp_path / "qc" / "dummy_qc.tsv"
    assert out.exists()
    lines = out.read_text().splitlines()
    assert lines[0] == "metric\tvalue"
    assert "one\t1" in lines[1:]
    assert "two\tdos" in lines[1:]
