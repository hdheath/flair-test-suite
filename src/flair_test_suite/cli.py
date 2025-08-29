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
from .validation import validate_stage_order
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
    """Yield Config objects assembled from a TSV manifest (base-first only).

    Required row format (tab-separated, paths resolved relative to the TSV):
        base_config.toml\tstage1.toml\tstage2.toml\t...

    The base config must define test_set_id.
    """
    with tsv_path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            if len(row) < 2:
                raise ValueError("Each TSV row must have at least a base config and one stage config")
            first = row[0].strip()
            rest = [c.strip() for c in row[1:] if c.strip()]

            def _resolve(p: str) -> Path:
                path = Path(p)
                if not path.is_absolute():
                    path = (tsv_path.parent / path).resolve()
                return path

            base_path = _resolve(first)
            if not base_path.exists() or base_path.suffix.lower() != ".toml":
                raise ValueError(
                    f"First column must be a TOML base config; got '{first}'"
                )
            stage_paths = [_resolve(p) for p in rest]
            if not stage_paths:
                raise ValueError(
                    "Each row must include at least one stage fragment after the base config"
                )
            cfg = load_config_fragments(base_path, stage_paths)
            # Remember source TSV for logging/diagnostics
            try:
                setattr(cfg, "_path", str(tsv_path))
            except Exception:
                pass
            rid = (
                getattr(cfg, "test_set_id", None)
                or getattr(cfg.run, "test_set_id", None)
                or getattr(cfg, "run_id", None)
                or getattr(cfg.run, "run_id", None)
            )
            if not rid:
                raise ValueError(
                    f"Base config '{base_path}' is missing test_set_id; add test_set_id under the top-level or [run] section."
                )
            yield cfg



def _load_cfg(cfg_path: Path):
    try:
        return load_config(cfg_path)
    except Exception as e:
        warnings.warn(f"Failed to load {cfg_path}: {e}", UserWarning)
        return None


def _prepare_logfile(cfg, run_id: str, work_dir: Path, seen_run_ids: set, quiet: bool) -> Path:
    # Logs are written under <work_dir>/<run_id>/run_summary.log
    log_dir = work_dir / run_id
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
            # Always display paths relative to the work_dir for readability
            out_path = format_path(pb.stage_dir, work_dir)
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
    inputs: Iterable[Union[Path, Config]]
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

        # Determine identifier from test_set_id
        run_id = (
            getattr(cfg, "test_set_id", None)
            or getattr(cfg.run, "test_set_id", None)
            or getattr(cfg, "run_id", None)
            or getattr(cfg.run, "run_id", None)
        )
        # Hardcode work_dir to ./outputs regardless of config
        work_dir = Path("./outputs")
        log_path = _prepare_logfile(cfg, run_id, work_dir, seen_run_ids, quiet=True)
        # Log which case index we are executing (1-based) for this TSV
        idx = getattr(cfg, "_case_idx", None)
        total = getattr(cfg, "_case_total", None)
        if idx and total:
            logger.info("Executing test case %s/%s from %s", idx, total, cfg_path)
        else:
            logger.info("Loading config: %s", cfg_path)

        # Validate TSV ordering and duplicates before attempting to run
        try:
            validate_stage_order(cfg)
        except Exception as e:
            click.echo(f"Invalid stage order: {e}", err=True)
            logger.error("Invalid stage order: %s", e)
            any_failed = True
            continue

        stage_order = topological_sort(cfg.run.stages)
        upstreams: Dict[str, PathBuilder] = {}
        # Print a single concise stdout message listing the stages to be executed
        click.echo(f"Run summary: {click.style(str(log_path), fg='cyan')}")
        try:
            stage_names = [st_cfg.name for st_cfg in stage_order]
        except Exception:
            stage_names = [str(s) for s in stage_order]
        print(f"Starting run '{run_id}' with {len(stage_order)} stage(s): {', '.join(stage_names)}")

        all_skipped = True
        for st_cfg in stage_order:
            failed, skipped = _execute_stage(
                st_cfg, cfg, run_id, work_dir, upstreams
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
def main(config_input: Path) -> None:
    """Entry point for the CLI."""
    if config_input.suffix.lower() in (".tsv", ".txt"):
        inputs = list(parse_cases_tsv(config_input))
        # Annotate each Config with its index to log in the run summary
        total = len(inputs)
        for i, cfg in enumerate(inputs, 1):
            try:
                setattr(cfg, "_case_idx", i)
                setattr(cfg, "_case_total", total)
            except Exception:
                pass
    else:
        raise click.UsageError("Only TSV manifests are supported. Provide a TSV with: base_config\tstage1\tstage2...")
    sys.exit(run_configs(inputs))


if __name__ == "__main__":  # pragma: no cover
    main()
