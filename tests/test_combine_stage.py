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
    # user supplied sample and manifest
    (data_dir / "s1.isoforms.bed").write_text("bed2")
    (data_dir / "s1.isoforms.fa").write_text("fa2")
    manifest_file = data_dir / "existing_manifest.tsv"
    manifest_file.write_text(
        f"1\tisoform\t{data_dir/'s1.isoforms.bed'}\t{data_dir/'s1.isoforms.fa'}\t\n"
    )

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
                    "flags": {"manifest": "existing_manifest.tsv"},
                }
            ],
        },
    }
    cfg = Cfg.model_validate(cfg_dict)

    collapse_pb = PathBuilder(repo_root / "outputs", "run1", "collapse", "sigC")
    collapse_dir = collapse_pb.stage_dir
    collapse_dir.mkdir(parents=True)
    collapse_bed = collapse_dir / "run1.isoforms.bed"
    collapse_bed.write_text("bedc")
    collapse_fa = collapse_dir / "run1.isoforms.fa"
    collapse_fa.write_text("fac")
    collapse_rm = collapse_dir / "run1.isoform.read.map.txt"
    collapse_rm.write_text("rmc")
    stage = CombineStage(
        cfg, run_id="run1", work_dir=repo_root / "outputs", upstreams={"collapse": collapse_pb}
    )

    cmd = stage.build_cmds()[0]
    assert cmd[:3] == ["flair", "combine", "--manifest"]
    assert cmd[3] == "manifest.tsv"
    assert cmd[4:6] == ["-o", "run1"]

    expected_hash = {collapse_bed, collapse_fa, collapse_rm, manifest_file}
    assert expected_hash.issubset(set(stage._hash_inputs))

    def fake_run_all(self, cmds, log_path, cwd):
        (cwd / "run1.bed").write_text("bed")
        (cwd / "run1.counts.tsv").write_text("counts")
        (cwd / "run1.fa").write_text("fa")
        (cwd / "run1.isoform.map.txt").write_text("map")
        return 0

    monkeypatch.setattr(StageBase, "_run_all", fake_run_all)

    pb = stage.run()
    manifest_path = pb.stage_dir / "manifest.tsv"
    assert manifest_path.exists()
    lines = manifest_path.read_text().splitlines()
    assert lines == [
        f"1\tisoform\t{data_dir/'s1.isoforms.bed'}\t{data_dir/'s1.isoforms.fa'}\t",
        f"2\tisoform\t{collapse_bed}\t{collapse_fa}\t{collapse_rm}",
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
    collapse_bed = collapse_dir / "run1.isoforms.bed"
    collapse_bed.write_text("bedc")

    stage = CombineStage(cfg, run_id="run1", work_dir=repo_root / "outputs", upstreams={"collapse": collapse_pb})

    cmd = stage.build_cmds()[0]
    assert cmd[:3] == ["flair", "combine", "--manifest"]
    assert cmd[3] == "manifest.tsv"
    assert cmd[4:6] == ["-o", "run1"]

    def fake_run_all(self, cmds, log_path, cwd):
        (cwd / "run1.bed").write_text("bed")
        (cwd / "run1.counts.tsv").write_text("counts")
        (cwd / "run1.fa").write_text("fa")
        (cwd / "run1.isoform.map.txt").write_text("map")
        return 0

    monkeypatch.setattr(StageBase, "_run_all", fake_run_all)

    pb = stage.run()
    manifest_path = pb.stage_dir / "manifest.tsv"
    assert manifest_path.exists()
    lines = manifest_path.read_text().splitlines()
    assert lines == [f"1\tisoform\t{collapse_bed}\t\t"]
