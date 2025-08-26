import sys
import types
from pathlib import Path
import json

# Make package importable
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.modules.pop("flair_test_suite", None)
sys.modules.pop("flair_test_suite.stages", None)

from flair_test_suite.stages import CollapseStage, TranscriptomeStage
from flair_test_suite.lib.reinstate import Reinstate


def _setup_stage(stage_cls, stage_name, tmp_path, monkeypatch, make_browser=True, region_tag=None):
    stage_dir = tmp_path / "out" / "run1" / stage_name / "sig"
    stage_dir.mkdir(parents=True)
    if region_tag:
        primary = stage_dir / f"{region_tag}.isoforms.bed"
    else:
        primary = stage_dir / "run1.isoforms.bed"
    primary.write_text("chr1\t0\t1\tread\t0\t+\n")

    def fake_collect(stage_dir, cfg, upstreams=None, out_dir=None):
        out = out_dir or stage_dir
        out.mkdir(parents=True, exist_ok=True)
        (out / "TED.tsv").write_text("col1\tcol2\n")
        if make_browser:
            br = out / "transcriptome_browser"
            br.mkdir()
            png = br / "dummy.png"
            png.write_text("png")
            (br / "region_map.json").write_text(json.dumps({"tag": str(png)}))

    import flair_test_suite.qc.ted as ted_mod
    monkeypatch.setattr(ted_mod, "collect", fake_collect)

    cfg = types.SimpleNamespace(run=types.SimpleNamespace(version="3.0.0"))
    stage = stage_cls(cfg, "run1", tmp_path, upstreams={})
    stage._run_qc(stage_dir, primary, runtime=None)
    (stage_dir / ".completed.json").write_text("{}")
    action = Reinstate.decide(stage_dir, primary, needs_qc=True, stage_name=stage_name)
    return action


def test_collapse_reinstate(tmp_path, monkeypatch):
    action = _setup_stage(CollapseStage, "collapse", tmp_path, monkeypatch)
    assert action == "skip"


def test_transcriptome_reinstate(tmp_path, monkeypatch):
    action = _setup_stage(TranscriptomeStage, "transcriptome", tmp_path, monkeypatch)
    assert action == "skip"


def test_missing_browser_triggers_qc(tmp_path, monkeypatch):
    action = _setup_stage(CollapseStage, "collapse", tmp_path, monkeypatch, make_browser=False, region_tag="chr1_0_10")
    assert action == "qc"

