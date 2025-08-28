import logging
import warnings
from pathlib import Path


def setup_run_logging(log_path: Path, mode: str = "a", quiet: bool = False) -> None:
    """Configure logging for a pipeline run.

    When ``quiet`` is true, only a file handler is configured so messages are
    written to ``log_path`` without appearing on the console.
    """
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    handlers: list[logging.Handler] = [logging.FileHandler(log_path, mode=mode)]
    if not quiet:
        handlers.append(logging.StreamHandler())

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
        handlers=handlers,
    )

    # Route Python warnings (warnings.warn) into the logging system so they
    # respect 'quiet' and land in the run summary log instead of stderr.
    logging.captureWarnings(True)
    # Optionally, downgrade noisy deprecation warnings
    try:
        warnings.filterwarnings("default")
    except Exception:
        pass
