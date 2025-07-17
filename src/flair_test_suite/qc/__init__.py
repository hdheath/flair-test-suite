"""
QC plug‑in registry.  Each stage registers its collector via @register().
"""

QC_REGISTRY: dict[str, callable] = {}

def register(stage_name: str):
    def _decorator(func):
        QC_REGISTRY[stage_name] = func
        return func
    return _decorator


# --- auto‑import built‑in collectors so they self‑register ------------
from . import align_qc  # noqa: F401  (<─ ensures decorator runs)
