import sys
from pathlib import Path

import pytest

# Make package importable
sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

# Ensure the real stages package is loaded (other tests may stub it)
import sys as _sys
_sys.modules.pop("flair_test_suite.stages", None)

from flair_test_suite.config_schema import Config as Cfg
from flair_test_suite.stages.combine import CombineStage
from flair_test_suite.stages.base import StageBase
from flair_test_suite.lib.paths import PathBuilder


@pytest.mark.skipif(not hasattr(StageBase, "run"), reason="StageBase missing run")
def test_combine_writes_manifest_and_builds_cmd(tmp_path, monkeypatch):
    repo_root = tmp_path
    data_dir = repo_root / "data"
    data_dir.mkdir()
    (data_dir / "s1.isoforms.bed").write_text("bed1")
    (data_dir / "s2.isoforms.bed").write_text("bed2")

    cfg_dict = {
        "run_id": "run1",
        "run": {
            "version": "1",
            "conda_env": "env",
            "work_dir": str(repo_root / "outputs"),
            "data_dir": str(data_dir),
            "stages": [
                {
                    "name": "combine",
                    "requires": ["collapse"],
                    "flags": {"manifest": ["s1.isoforms.bed", "s2.isoforms.bed"]},
                }
            ],
        },
    }
    cfg = Cfg.model_validate(cfg_dict)

    collapse_pb = PathBuilder(repo_root / "outputs", "run1", "collapse", "sigC")
    collapse_dir = collapse_pb.stage_dir
    collapse_dir.mkdir(parents=True)
    collapse_out = collapse_dir / "run1.isoforms.bed"
    collapse_out.write_text("bedc")
    stage = CombineStage(cfg, run_id="run1", work_dir=repo_root / "outputs", upstreams={"collapse": collapse_pb})

    cmd = stage.build_cmd()
    assert cmd[:3] == ["flair", "combine", "-i"]
    assert cmd[3] == "manifest.tsv"
    assert cmd[4:6] == ["-o", "run1"]

    def fake_run_all(self, cmds, log_path, cwd):
        (cwd / "run1_combined.isoforms.bed").write_text("bed")
        (cwd / "run1_combined.isoforms.gtf").write_text("gtf")
        return 0

    monkeypatch.setattr(StageBase, "_run_all", fake_run_all)

    pb = stage.run()
    manifest_path = pb.stage_dir / "manifest.tsv"
    assert manifest_path.exists()
    lines = manifest_path.read_text().strip().splitlines()
    assert lines == [
        str(collapse_out),
        str(data_dir / "s1.isoforms.bed"),
        str(data_dir / "s2.isoforms.bed"),
    ]


@pytest.mark.skipif(not hasattr(StageBase, "run"), reason="StageBase missing run")
def test_combine_defaults_to_collapse_output(tmp_path, monkeypatch):
    repo_root = tmp_path
    data_dir = repo_root / "data"
    data_dir.mkdir()

    cfg_dict = {
        "run_id": "run1",
        "run": {
            "version": "1",
            "conda_env": "env",
            "work_dir": str(repo_root / "outputs"),
            "data_dir": str(data_dir),
            "stages": [
                {
                    "name": "combine",
                    "requires": ["collapse"],
                }
            ],
        },
    }
    cfg = Cfg.model_validate(cfg_dict)

    collapse_pb = PathBuilder(repo_root / "outputs", "run1", "collapse", "sigC")
    collapse_dir = collapse_pb.stage_dir
    collapse_dir.mkdir(parents=True)
    collapse_out = collapse_dir / "run1.isoforms.bed"
    collapse_out.write_text("bedc")

    stage = CombineStage(cfg, run_id="run1", work_dir=repo_root / "outputs", upstreams={"collapse": collapse_pb})

    cmd = stage.build_cmd()
    assert cmd[:3] == ["flair", "combine", "-i"]
    assert cmd[3] == "manifest.tsv"
    assert cmd[4:6] == ["-o", "run1"]

    def fake_run_all(self, cmds, log_path, cwd):
        (cwd / "run1_combined.isoforms.bed").write_text("bed")
        (cwd / "run1_combined.isoforms.gtf").write_text("gtf")
        return 0

    monkeypatch.setattr(StageBase, "_run_all", fake_run_all)

    pb = stage.run()
    manifest_path = pb.stage_dir / "manifest.tsv"
    assert manifest_path.exists()
    lines = manifest_path.read_text().strip().splitlines()
    assert lines == [str(collapse_out)]
