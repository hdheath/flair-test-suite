from pathlib import Path
from click.testing import CliRunner
import sys
import types
import tomllib

sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

config_loader = types.ModuleType("flair_test_suite.config_loader")


class StageConfig:
    def __init__(self, name, requires=None):
        self.name = name
        self.requires = requires or []


class RunConfig:
    def __init__(self, stages, **kwargs):
        self.stages = stages
        for k, v in kwargs.items():
            setattr(self, k, v)


class Config:
    def __init__(self, run_id, run):
        self.run_id = run_id
        self.run = run


def load_config(path):
    with open(path, "rb") as fh:
        data = tomllib.load(fh)
    run_data = data["run"]
    stages = [StageConfig(**s) for s in run_data["stages"]]
    run = RunConfig(stages=stages, **{k: v for k, v in run_data.items() if k != "stages"})
    return Config(run_id=data.get("run_id"), run=run)


config_loader.load_config = load_config
sys.modules["flair_test_suite.config_loader"] = config_loader

stages_module = types.ModuleType("flair_test_suite.stages")
stages_module.STAGE_REGISTRY = {}
sys.modules["flair_test_suite.stages"] = stages_module

from flair_test_suite import cli
from flair_test_suite.lib import PathBuilder


class DummyStage:
    name = "dummy"

    def __init__(self, cfg, run_id, work_dir, upstreams):
        self.cfg = cfg
        self.run_id = run_id
        self.work_dir = Path(work_dir)
        self.action = "run"

    def build_cmd(self):
        return None

    def run(self):
        pb = PathBuilder(self.work_dir, self.run_id, self.name, "sig")
        pb.stage_dir.mkdir(parents=True, exist_ok=True)
        return pb


def _write_config(tmp_path: Path) -> Path:
    """Create a minimal TOML config pointing to DummyStage."""
    work_dir = tmp_path / "work"
    cfg_text = f"""
run_id = "demo"
[run]
version = "1"
conda_env = "env"
work_dir = "{work_dir}"
data_dir = "{tmp_path}"
stages = [{{name = "dummy"}}]
"""
    cfg_path = tmp_path / "config.toml"
    cfg_path.write_text(cfg_text)
    return cfg_path


def test_iter_config_paths(tmp_path: Path):
    cfg1 = tmp_path / "a.toml"
    cfg2 = tmp_path / "sub" / "b.toml"
    cfg1.touch()
    cfg2.parent.mkdir()
    cfg2.touch()
    tsv = tmp_path / "list.tsv"
    tsv.write_text(f"{cfg1.name}\n# comment\nsub/b.toml\n")

    paths = list(cli._iter_config_paths(tsv))
    assert paths == [cfg1, cfg2]


def test_run_configs_executes_stage(tmp_path: Path, monkeypatch):
    cfg = _write_config(tmp_path)
    monkeypatch.setitem(cli.STAGE_REGISTRY, "dummy", DummyStage)
    code = cli.run_configs([cfg])
    assert code == 0
    assert (tmp_path / "work" / "demo" / "dummy" / "sig").exists()


def test_cli_main(tmp_path: Path, monkeypatch):
    cfg = _write_config(tmp_path)
    monkeypatch.setitem(cli.STAGE_REGISTRY, "dummy", DummyStage)
    runner = CliRunner()
    result = runner.invoke(cli.main, [str(cfg)])
    assert result.exit_code == 0
    assert (tmp_path / "work" / "demo" / "dummy" / "sig").exists()
