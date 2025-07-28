import logging
from pathlib import Path

def setup_run_logging(log_path: Path, mode: str = "a"):
    """Configure logging for a pipeline run, overwriting or appending as needed."""
    # Remove all existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    # Set up new handlers
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_path, mode=mode)
        ]
    )