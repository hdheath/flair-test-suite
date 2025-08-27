# src/flair_test_suite/cli.py

import sys
import warnings
import logging
import csv
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, Tuple, Union

import click

from .config_loader import load_config, load_config_fragments
from .config_schema import Config
from .stages import STAGE_REGISTRY
from .lib import topological_sort
from .lib.logging import setup_run_logging
from .lib.paths import PathBuilder, format_path


logger = logging.getLogger(__name__)


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


def parse_cases_tsv(tsv_path: Path) -> Iterator[Config]:
    """Yield Config objects assembled from a TSV manifest.

    Each non-empty line must contain at least a run ID and base config path
    followed by zero or more stage config paths. Paths are resolved relative
    to the TSV location.
    """
    with tsv_path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            if len(row) < 2:
                raise ValueError("Each row must have run_id and base config path")
            run_id = row[0].strip()
            base_rel = row[1].strip()
            stage_rel = [c.strip() for c in row[2:] if c.strip()]

            def _resolve(p: str) -> Path:
                path = Path(p)
                if not path.is_absolute():
                    path = (tsv_path.parent / path).resolve()
                return path

            base_path = _resolve(base_rel)
            stage_paths = [_resolve(p) for p in stage_rel]

            cfg = load_config_fragments(base_path, stage_paths)
            cfg.run_id = run_id
            yield cfg



def _load_cfg(cfg_path: Path):
    try:
        return load_config(cfg_path)
    except Exception as e:
        warnings.warn(f"Failed to load {cfg_path}: {e}", UserWarning)
        return None


def _prepare_logfile(cfg, run_id: str, work_dir: Path, seen_run_ids: set, quiet: bool) -> Path:
    case_name = getattr(cfg, "case_name", None)
    log_dir = work_dir
    if case_name:
        log_dir = log_dir / case_name
    log_dir = log_dir / run_id
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "run_summary.log"
    mode = "a" if run_id in seen_run_ids else "w"
    setup_run_logging(log_path, mode=mode, quiet=quiet)
    seen_run_ids.add(run_id)
    logger.info("Run summary: %s", format_path(log_path, work_dir))
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

    logger.info("Starting %s", st_cfg.name)
    click.echo(f"Running {st_cfg.name}")

    try:
        pb = stage_instance.run()
    except Exception as e:
        logger.error("Stage %s failed: %s", st_cfg.name, e)
        click.echo(f"Stage {st_cfg.name} failed: {e}", err=True)
        return True, True  # failed, skipped
    else:
        upstreams[st_cfg.name] = pb
        skipped = getattr(stage_instance, "action", None) == "skip"
        if pb is not None:
            path_root = None if absolute_paths else work_dir
            out_path = format_path(pb.stage_dir, path_root)
            logger.info("\u2713 Done %s – outputs: %s", st_cfg.name, out_path)
        else:
            logger.warning("Stage %s produced no outputs.", st_cfg.name)
        if skipped:
            click.echo(f"Skipping {st_cfg.name}")
            logger.info("Skipped %s", st_cfg.name)
        else:
            click.echo(f"Completed {st_cfg.name}")
            logger.info("Completed %s", st_cfg.name)
        return False, skipped


# ────────────────────────── main orchestration ──────────────────────────
def run_configs(
    inputs: Iterable[Union[Path, Config]], absolute_paths: bool = False
) -> int:
    """Run one or more configuration files and return an exit code."""

    any_ran = False
    any_failed = False
    all_configs_skipped = True
    seen_run_ids: set[str] = set()

    for item in inputs:
        if isinstance(item, Path):
            cfg_path = item
            if not cfg_path.exists():
                warnings.warn(f"Config file not found: {cfg_path}", UserWarning)
                continue
            cfg = _load_cfg(cfg_path)
            if not cfg:
                continue
        else:
            cfg = item
            cfg_path = getattr(cfg, "_path", "<tsv>")

        run_id = getattr(cfg, "run_id", None) or getattr(cfg.run, "run_id", None)
        work_dir = Path(cfg.run.work_dir)
        log_path = _prepare_logfile(cfg, run_id, work_dir, seen_run_ids, quiet=True)
        logger.info("Loading config: %s", cfg_path)

        stage_order = topological_sort(cfg.run.stages)
        upstreams: Dict[str, PathBuilder] = {}
        click.echo(f"Starting run '{run_id}' with {len(stage_order)} stage(s)")
        click.echo(f"Run summary: {click.style(str(log_path), fg='cyan')}")
        # Also print a concise stdout message listing the stages to be executed
        try:
            stage_names = [st_cfg.name for st_cfg in stage_order]
        except Exception:
            stage_names = [str(s) for s in stage_order]
        print(f"Starting run '{run_id}' with {len(stage_order)} stage(s): {', '.join(stage_names)}")

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
        logger.info("All configurations are up to date; nothing to do.")
        return 0
    if any_failed:
        logger.error(
            "Pipeline failed: at least one stage did not complete successfully."
        )
        return 1
    if any_ran:
        logger.info("All stages completed successfully.")
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
    if config_input.suffix.lower() in (".tsv", ".txt"):
        lines = [
            l
            for l in config_input.read_text().splitlines()
            if l.strip() and not l.startswith("#")
        ]
        if lines and lines[0].count("\t") >= 2:
            inputs = list(parse_cases_tsv(config_input))
        else:
            inputs = list(_iter_config_paths(config_input))
    else:
        inputs = [config_input]
    sys.exit(run_configs(inputs, absolute_paths))


if __name__ == "__main__":  # pragma: no cover
    main()
