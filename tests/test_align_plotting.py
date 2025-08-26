from pathlib import Path
import json
import types
import sys

import pytest

sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

pytest.importorskip("pysam")

from flair_test_suite.plotting.align_histograms import generate_histograms
from flair_test_suite.qc import align_qc


def test_generate_histograms(tmp_path):
    out = tmp_path / "plots"
    mapping = generate_histograms([1, 2], [50.0, 60.0], [100, 200], out)
    assert (out / mapping["mapq"]).exists()
    assert (out / mapping["identity"]).exists()
    assert (out / mapping["length"]).exists()


def test_align_qc_delegates_plotting(tmp_path, monkeypatch):
    bam = tmp_path / "in.bam"
    bam.touch()
    bed = tmp_path / "in.bed"
    bed.write_text("chr1\t0\t1\t.\n")

    # stub external calls
    def fake_run(args, check=True, stdout=None, stderr=None):
        if "-o" in args:
            out = Path(args[args.index("-o") + 1])
            out.write_text("MAPQ\t10\t1\nRL\t100\t1\n")
        return None

    monkeypatch.setattr(
        align_qc,
        "subprocess",
        types.SimpleNamespace(run=fake_run, CalledProcessError=Exception, DEVNULL=None),
    )
    monkeypatch.setattr(align_qc, "iter_primary", lambda *a, **k: [])
    monkeypatch.setattr(align_qc, "count_unique_junctions", lambda *a, **k: 0)
    monkeypatch.setattr(align_qc, "count_splice_junction_motifs", lambda *a, **k: {})
    class DummyAF:
        def __enter__(self):
            return object()
        def __exit__(self, *a):
            return False

    monkeypatch.setattr(align_qc, "pysam", types.SimpleNamespace(AlignmentFile=lambda *a, **k: DummyAF()))

    called = {}

    def fake_generate(mapq, ident, length, out_dir):
        called["out_dir"] = out_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        return {"mapq": "a.png", "identity": "b.png", "length": "c.png"}

    monkeypatch.setattr(align_qc, "generate_histograms", fake_generate)

    metrics = align_qc.collect(bam, tmp_path, n_input_reads=10, genome_fa="gen.fa", runtime_sec=1.0, read_count_method="estimated")
    qc_dir = tmp_path / "qc"
    assert called["out_dir"] == qc_dir
    manifest = qc_dir / "align_plot_manifest.json"
    assert manifest.exists()
    data = json.loads(manifest.read_text())
    assert data["mapq"] == "a.png"
    assert metrics["read_count_method"] == "estimated"
