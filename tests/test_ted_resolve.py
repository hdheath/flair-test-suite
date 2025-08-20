import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1] / 'src'))

from flair_test_suite.qc.ted import _resolve


def test_resolve_expands_env_and_user(tmp_path, monkeypatch):
    monkeypatch.setenv('DATA_HOME', str(tmp_path))
    # env var expansion
    p = _resolve('${DATA_HOME}/file.bed', Path('/unused'))
    assert p == tmp_path / 'file.bed'

    # user expansion
    monkeypatch.setenv('HOME', str(tmp_path))
    p2 = _resolve('~/another.bed', Path('/unused'))
    assert p2 == tmp_path / 'another.bed'
