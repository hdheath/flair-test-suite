from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

from flair_test_suite.cli import parse_cases_tsv


def test_parse_cases_tsv(tmp_path: Path):
    tsv = tmp_path / "cases.tsv"
    tsv.write_text(
        """key\tvalue
test_set_id\trun1
version\t1
flair_env\tenv
data_dir\tdata
reads_file\treads.fa

stage\tflags
align\t--foo, --bar 1
collapse\t
"""
    )

    cfgs = list(parse_cases_tsv(tsv))
    assert len(cfgs) == 1
    cfg = cfgs[0]
    assert cfg.test_set_id == "run1"
    assert [s.name for s in cfg.run.stages] == ["align", "collapse"]
