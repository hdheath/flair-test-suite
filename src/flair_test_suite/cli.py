# src/flair_test_suite/cli.py

import sys
import warnings
import logging
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, Tuple

import click

from .config_loader import load_config
from .stages import STAGE_REGISTRY
from .lib import PathBuilder, topological_sort
from .lib.logging import setup_run_logging


# ────────────────────────── helpers ──────────────────────────
def _iter_config_paths(tsv_path: Path) -> Iterator[Path]:
    """Yield absolute Paths to config files listed in the TSV."""
    for raw in tsv_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        rel = line.split("\t")[0].strip()
        cfg_path = Path(rel)
        if not cfg_path.is_absolute():
            cfg_path = (tsv_path.parent / cfg_path).resolve()
        yield cfg_path


def _setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.StreamHandler()],
    )


def _load_cfg(cfg_path: Path):
    try:
        return load_config(cfg_path)
    except Exception as e:
        warnings.warn(f"Failed to load {cfg_path}: {e}", UserWarning)
        return None


def _prepare_logfile(cfg, run_id: str, work_dir: Path, seen_run_ids: set) -> None:
    log_path = work_dir / run_id / "run_summary.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if run_id not in seen_run_ids else "a"
    seen_run_ids.add(run_id)
    setup_run_logging(log_path, mode)


def _execute_stage(st_cfg, cfg, run_id: str, work_dir: Path, upstreams: Dict[str, PathBuilder]) -> Tuple[bool, bool]:
    """Execute one stage. Returns (failed, skipped)."""
    StageCls = STAGE_REGISTRY[st_cfg.name]
    stage_instance = StageCls(cfg, run_id, work_dir, upstreams)
    stage_instance.build_cmd()

    try:
        pb = stage_instance.run()
    except Exception as e:
        logging.error(f"Stage {st_cfg.name} failed: {e}")
        return True, True   # failed, skipped
    else:
        upstreams[st_cfg.name] = pb
        if pb is not None:
            logging.info(f"✓ Done {st_cfg.name} – outputs: {pb.stage_dir}")
        else:
            logging.warning(f"Stage {st_cfg.name} produced no outputs.")
        skipped = getattr(stage_instance, "action", None) == "skip"
        return False, skipped


# ────────────────────────── main orchestration ──────────────────────────
def run_configs(inputs: Iterable[Path]) -> int:
    """Run one or more configuration files and return an exit code."""
    _setup_logging()

    any_ran = False
    any_failed = False
    all_configs_skipped = True
    seen_run_ids: set[str] = set()

    for cfg_path in inputs:
        if not cfg_path.exists():
            warnings.warn(f"Config file not found: {cfg_path}", UserWarning)
            continue

        logging.info(f"Loading config: {cfg_path}")
        cfg = _load_cfg(cfg_path)
        if not cfg:
            continue

        run_id = getattr(cfg, "run_id", None) or getattr(cfg.run, "run_id", None)
        work_dir = Path(cfg.run.work_dir)
        _prepare_logfile(cfg, run_id, work_dir, seen_run_ids)

        stage_order = topological_sort(cfg.run.stages)
        upstreams: Dict[str, PathBuilder] = {}
        logging.info(f"Executing run '{run_id}' with {len(stage_order)} stage(s)")

        all_skipped = True
        for st_cfg in stage_order:
            failed, skipped = _execute_stage(st_cfg, cfg, run_id, work_dir, upstreams)
            if failed:
                any_failed = True
                break
            if not skipped:
                all_skipped = False
                any_ran = True

        if not all_skipped:
            all_configs_skipped = False

    # summary & exit code
    if not seen_run_ids:
        raise click.UsageError("No valid configuration files were executed.")

    if all_configs_skipped:
        logging.info("All configurations are up to date; nothing to do.")
        return 0
    if any_failed:
        logging.error("Pipeline failed: at least one stage did not complete successfully.")
        return 1
    if any_ran:
        logging.info("All stages completed successfully.")
    return 0


@click.command()
@click.argument("config_input", type=click.Path(exists=True, path_type=Path))
def main(config_input: Path) -> None:
    """Entry point for the CLI."""
    inputs = _iter_config_paths(config_input) if config_input.suffix.lower() in (".tsv", ".txt") else [config_input]
    sys.exit(run_configs(inputs))


if __name__ == "__main__":  # pragma: no cover
    main()
