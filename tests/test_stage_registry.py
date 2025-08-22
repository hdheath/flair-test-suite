from pathlib import Path
import types
import sys
import importlib
import pytest

# Ensure package imports resolve without installation
sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

# Skip if heavy QC deps needed by stage imports are missing
pytest.importorskip("pysam")
pytest.importorskip("matplotlib")

# Reload the real stages package in case other tests stubbed it
sys.modules.pop("flair_test_suite.stages", None)
stages_pkg = importlib.import_module("flair_test_suite.stages")
from flair_test_suite.stages.base import StageBase
STAGE_REGISTRY = stages_pkg.STAGE_REGISTRY


def test_stage_registry_contains_expected_classes():
    """Ensure shipped stages are registered with their names."""
    expected = {"align", "correct", "regionalize", "collapse", "combine", "transcriptome", "quantify"}
    assert expected <= STAGE_REGISTRY.keys()

    dummy_cfg = types.SimpleNamespace(run=types.SimpleNamespace(data_dir=".", conda_env="env"))
    for name in expected:
        cls = STAGE_REGISTRY[name]
        assert issubclass(cls, StageBase)
        stage = cls(dummy_cfg, run_id="demo", work_dir=Path("/tmp"), upstreams={})
        assert stage.name == name
        assert hasattr(stage, "primary_output_key")
