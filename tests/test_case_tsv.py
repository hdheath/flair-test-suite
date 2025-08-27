from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

from flair_test_suite.cli import parse_cases_tsv


def test_parse_cases_tsv(tmp_path: Path):
    base = tmp_path / "base.toml"
    base.write_text(
        """
case_name = "demo"
[run]
version = "1"
conda_env = "env"
work_dir = "work"
data_dir = "data"
"""
    )
    align = tmp_path / "align.toml"
    align.write_text(
        """
[[run.stages]]
name = "align"
"""
    )
    collapse = tmp_path / "collapse.toml"
    collapse.write_text(
        """
[[run.stages]]
name = "collapse"
requires = ["align"]
"""
    )
    tsv = tmp_path / "cases.tsv"
    tsv.write_text(
        f"run1\t{base.name}\t{align.name}\t{collapse.name}\n"
    )

    cfgs = list(parse_cases_tsv(tsv))
    assert len(cfgs) == 1
    cfg = cfgs[0]
    assert cfg.run_id == "run1"
    assert cfg.case_name == "demo"
    assert [s.name for s in cfg.run.stages] == ["align", "collapse"]

