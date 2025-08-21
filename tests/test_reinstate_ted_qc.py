import sys
import types
from pathlib import Path

# Make package importable
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.modules.pop("flair_test_suite", None)
sys.modules.pop("flair_test_suite.stages", None)

from flair_test_suite.stages import CollapseStage, TranscriptomeStage
from flair_test_suite.lib.reinstate import Reinstate


def _setup_stage(stage_cls, stage_name, tmp_path, monkeypatch, make_browser=True):
    stage_dir = tmp_path / "out" / "run1" / stage_name / "sig"
    stage_dir.mkdir(parents=True)
    primary = stage_dir / "run1.isoforms.bed"
    primary.write_text("chr1\t0\t1\tread\t0\t+\n")

    def fake_collect(stage_dir, cfg, upstreams=None):
        (stage_dir / "TED.tsv").write_text("col1\tcol2\n")
        if make_browser:
            br = stage_dir / "transcriptome_browser"
            br.mkdir()
            (br / "dummy.png").write_text("png")

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
    action = _setup_stage(CollapseStage, "collapse", tmp_path, monkeypatch, make_browser=False)
    assert action == "qc"

