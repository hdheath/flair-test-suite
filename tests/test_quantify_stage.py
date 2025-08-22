from pathlib import Path
import sys
import pytest

# Ensure package importable
sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

# Reset stage registry
import sys as _sys
_sys.modules.pop("flair_test_suite.stages", None)

from flair_test_suite.config_schema import Config as Cfg
from flair_test_suite.stages.quantify import QuantifyStage
from flair_test_suite.stages.base import StageBase
from flair_test_suite.lib.paths import PathBuilder


@pytest.mark.skipif(not hasattr(StageBase, "run"), reason="StageBase missing run")
def test_quantify_uses_manifest_and_cmd(tmp_path, monkeypatch):
    repo_root = tmp_path
    data_dir = repo_root / "data"
    data_dir.mkdir()
    reads = data_dir / "reads.fq"
    reads.write_text("@r1\nAAAA\n+\n!!!!\n")
    manifest = data_dir / "reads_manifest.tsv"
    manifest.write_text(f"run1\tcondition1\tbatch1\t{reads}\n")

    cfg_dict = {
        "run_id": "run1",
        "run": {
            "version": "1",
            "conda_env": "env",
            "work_dir": str(repo_root / "outputs"),
            "data_dir": str(data_dir),
            "stages": [
                {
                    "name": "quantify",
                    "requires": ["collapse"],
                    "flags": {"manifest": "reads_manifest.tsv"},
                }
            ],
        },
    }
    cfg = Cfg.model_validate(cfg_dict)

    collapse_pb = PathBuilder(repo_root / "outputs", "run1", "collapse", "sigC")
    collapse_dir = collapse_pb.stage_dir
    collapse_dir.mkdir(parents=True)
    iso_fa = collapse_dir / "run1.isoforms.fa"
    iso_fa.write_text(">iso1\nAAAA\n")

    stage = QuantifyStage(
        cfg, run_id="run1", work_dir=repo_root / "outputs", upstreams={"collapse": collapse_pb}
    )

    cmd = stage.build_cmd()
    assert cmd[:3] == ["flair", "quantify", "-r"]
    assert cmd[3] == "reads_manifest.tsv"
    assert "-i" in cmd and str(iso_fa) in cmd

    expected_hash = {iso_fa, reads, manifest}
    assert expected_hash.issubset(set(stage._hash_inputs))

    def fake_run_all(self, cmds, log_path, cwd):
        (cwd / "run1.isoform.tpm.tsv").write_text("isoform\ts1\niso1\t1\n")
        (cwd / "run1.gene.counts.tsv").write_text("gene\ts1\nG1\t1\n")
        return 0

    monkeypatch.setattr(StageBase, "_run_all", fake_run_all)

    pb = stage.run()
    manifest_path = pb.stage_dir / "reads_manifest.tsv"
    assert manifest_path.exists()
    lines = manifest_path.read_text().splitlines()
    assert lines == [f"run1\tcondition1\tbatch1\t{reads}"]


@pytest.mark.skipif(not hasattr(StageBase, "run"), reason="StageBase missing run")
def test_quantify_requires_manifest(tmp_path):
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
            "stages": [{"name": "quantify", "requires": ["collapse"]}],
        },
    }
    cfg = Cfg.model_validate(cfg_dict)

    collapse_pb = PathBuilder(repo_root / "outputs", "run1", "collapse", "sigC")
    collapse_dir = collapse_pb.stage_dir
    collapse_dir.mkdir(parents=True)
    (collapse_dir / "run1.isoforms.fa").write_text(">iso1\nAAAA\n")

    stage = QuantifyStage(
        cfg, run_id="run1", work_dir=repo_root / "outputs", upstreams={"collapse": collapse_pb}
    )

    with pytest.raises(RuntimeError):
        stage.run()
