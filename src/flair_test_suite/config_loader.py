import tomli
from pathlib import Path
from .config_schema import Config

def load_config(path: str | Path) -> Config:
    with open(path, "rb") as fh:
        data = tomli.load(fh)
    return Config.parse_obj(data)