# src/flair_test_suite/cli.py

import sys
import warnings
import logging
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, Tuple
from datetime import datetime
import shutil

import click

from .config_loader import load_config
from .stages import STAGE_REGISTRY
from .lib import topological_sort
from .lib.logging import setup_run_logging
from .lib.paths import PathBuilder, format_path


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



def _load_cfg(cfg_path: Path):
    try:
        return load_config(cfg_path)
    except Exception as e:
        warnings.warn(f"Failed to load {cfg_path}: {e}", UserWarning)
        return None


def _prepare_logfile(cfg, run_id: str, work_dir: Path, seen_run_ids: set, quiet: bool) -> Path:
    log_dir = work_dir / run_id
    log_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    log_path = log_dir / f"run_summary_{ts}.log"
    setup_run_logging(log_path, mode="w", quiet=quiet)
    latest = log_dir / "run_summary.log"
    try:
        if latest.exists() or latest.is_symlink():
            latest.unlink()
        latest.symlink_to(log_path.name)
    except OSError:
        shutil.copy2(log_path, latest)
    seen_run_ids.add(run_id)
    logging.info(f"Run summary: {format_path(log_path, work_dir)}")
    return log_path


def _execute_stage(
    st_cfg,
    cfg,
    run_id: str,
    work_dir: Path,
    upstreams: Dict[str, PathBuilder],
    absolute_paths: bool,
) -> Tuple[bool, bool]:
    """Execute one stage. Returns (failed, skipped)."""
    StageCls = STAGE_REGISTRY[st_cfg.name]
    stage_instance = StageCls(cfg, run_id, work_dir, upstreams)
    try:
        stage_instance._build_and_normalize_cmds()
    except Exception:
        pass

    logging.info(f"Starting {st_cfg.name}")
    click.echo(f"Running {st_cfg.name}")

    try:
        pb = stage_instance.run()
    except Exception as e:
        logging.error(f"Stage {st_cfg.name} failed: {e}")
        click.echo(f"Stage {st_cfg.name} failed: {e}", err=True)
        return True, True  # failed, skipped
    else:
        upstreams[st_cfg.name] = pb
        skipped = getattr(stage_instance, "action", None) == "skip"
        if pb is not None:
            path_root = None if absolute_paths else work_dir
            out_path = format_path(pb.stage_dir, path_root)
            logging.info(f"✓ Done {st_cfg.name} – outputs: {out_path}")
        else:
            logging.warning(f"Stage {st_cfg.name} produced no outputs.")
        if skipped:
            click.echo(f"Skipping {st_cfg.name}")
        else:
            click.echo(f"Completed {st_cfg.name}")
        return False, skipped


# ────────────────────────── main orchestration ──────────────────────────
def run_configs(inputs: Iterable[Path], absolute_paths: bool = False) -> int:
    """Run one or more configuration files and return an exit code."""

    any_ran = False
    any_failed = False
    all_configs_skipped = True
    seen_run_ids: set[str] = set()

    for cfg_path in inputs:
        if not cfg_path.exists():
            warnings.warn(f"Config file not found: {cfg_path}", UserWarning)
            continue

        cfg = _load_cfg(cfg_path)
        if not cfg:
            continue

        run_id = getattr(cfg, "run_id", None) or getattr(cfg.run, "run_id", None)
        work_dir = Path(cfg.run.work_dir)
        log_path = _prepare_logfile(cfg, run_id, work_dir, seen_run_ids, quiet=True)
        logging.info(f"Loading config: {cfg_path}")

        stage_order = topological_sort(cfg.run.stages)
        upstreams: Dict[str, PathBuilder] = {}
        click.echo(f"Starting run '{run_id}' with {len(stage_order)} stage(s)")
        click.echo(f"Run summary: {click.style(str(log_path), fg='cyan')}")

        all_skipped = True
        for st_cfg in stage_order:
            failed, skipped = _execute_stage(
                st_cfg, cfg, run_id, work_dir, upstreams, absolute_paths
            )
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
        logging.error(
            "Pipeline failed: at least one stage did not complete successfully."
        )
        return 1
    if any_ran:
        logging.info("All stages completed successfully.")
    return 0


@click.command()
@click.argument("config_input", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--absolute-paths/--relative-paths",
    default=False,
    help="Show absolute paths in output",
)
def main(config_input: Path, absolute_paths: bool) -> None:
    """Entry point for the CLI."""
    inputs = (
        _iter_config_paths(config_input)
        if config_input.suffix.lower() in (".tsv", ".txt")
        else [config_input]
    )
    sys.exit(run_configs(inputs, absolute_paths))


if __name__ == "__main__":  # pragma: no cover
    main()
