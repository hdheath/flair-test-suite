## src/config_loader.py
import os

try:
    import tomllib
except ImportError:
    import toml as tomllib


def load_config(path: str) -> dict:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Config file not found: {path}")
    with open(path, "rb") as f:
        return tomllib.load(f)


def validate_config(cfg: dict):
    required = ["run_id", "input_root", "data_dir", "run"]
    missing = [k for k in required if k not in cfg]
    if missing:
        raise KeyError(f"Missing config keys: {missing}")


class Config:
    def __init__(self, raw: dict):
        validate_config(raw)
        self.raw        = raw
        self.run_id     = raw["run_id"]
        self.summary    = raw.get("summary", "")
        self.input_root = raw["input_root"]
        self.data       = raw["data_dir"]
        self.run        = raw["run"]

    @property
    def data_dir(self) -> str:
        return os.path.join(self.input_root, self.data["name"])