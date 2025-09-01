from pathlib import Path
from click.testing import CliRunner
import sys

sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

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


def _write_cases(tmp_path: Path) -> Path:
    """Create a minimal two-section TSV with a single dummy stage."""
    cases = tmp_path / "cases.tsv"
    cases.write_text(
        f"key\tvalue\n"
        f"test_set_id\tdemo\n"
        f"version\t1\n"
        f"flair_env\tenv\n"
        f"data_dir\t{tmp_path}\n"
        f"reads_file\treads.fa\n\n"
        f"stage\tflags\n"
        f"dummy\t\n"
    )
    return cases


def test_iter_config_paths(tmp_path: Path):
    # Still works for path lists (legacy helper, used elsewhere)
    p1 = tmp_path / "a.toml"; p1.touch()
    p2 = tmp_path / "sub" / "b.toml"; p2.parent.mkdir(); p2.touch()
    tsv = tmp_path / "list.tsv"; tsv.write_text(f"{p1.name}\n# comment\nsub/b.toml\n")
    paths = list(cli._iter_config_paths(tsv))
    assert paths == [p1, p2]


def test_run_configs_executes_stage(tmp_path: Path, monkeypatch):
    cases = _write_cases(tmp_path)
    monkeypatch.setitem(cli.STAGE_REGISTRY, "dummy", DummyStage)
    cfgs = list(cli.parse_cases_tsv(cases))
    code = cli.run_configs(cfgs)
    assert code == 0
    assert (Path("outputs") / "demo" / "dummy" / "sig").exists()


def test_cli_main(tmp_path: Path, monkeypatch):
    cases = _write_cases(tmp_path)
    monkeypatch.setitem(cli.STAGE_REGISTRY, "dummy", DummyStage)
    runner = CliRunner()
    result = runner.invoke(cli.main, [str(cases)])
    assert result.exit_code == 0
    assert (Path("outputs") / "demo" / "dummy" / "sig").exists()
