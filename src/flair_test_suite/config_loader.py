#src/flair_test_suite/config.py
import toml
from pathlib import Path
from types import SimpleNamespace

def _to_ns(obj):
    """
    Recursively convert dicts → SimpleNamespace and lists of dicts → lists of namespaces.
    """
    if isinstance(obj, dict):
        return SimpleNamespace(**{k: _to_ns(v) for k, v in obj.items()})
    elif isinstance(obj, list):
        return [_to_ns(v) for v in obj]
    else:
        return obj

def load_config(path):
    """
    Load a TOML file from `path` and return a namespaced config object.
    Usage:
        cfg = load_config("config/config_template.toml")
        print(cfg.run.version)
        print(cfg.run.steps.align.flags.threads)
    """
    text = Path(path).read_text(encoding="utf-8")
    raw  = toml.loads(text)
    return _to_ns(raw)
