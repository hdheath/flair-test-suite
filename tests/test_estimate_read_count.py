from pathlib import Path
import gzip
import sys

sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

from flair_test_suite.stages.stage_utils import estimate_read_count

def _write_fastq(path: Path, n: int):
    with gzip.open(path, "wt") as fh:
        for i in range(n):
            fh.write(f"@r{i}\nAAA\n+\nIII\n")


def test_estimate_read_count_exact(tmp_path):
    fq = tmp_path / "reads.fastq.gz"
    _write_fastq(fq, 10)
    count, exact = estimate_read_count(fq, max_reads=50)
    assert count == 10 and exact


def test_estimate_read_count_estimated(tmp_path):
    fq = tmp_path / "reads.fastq.gz"
    _write_fastq(fq, 1000)
    count, exact = estimate_read_count(fq, max_reads=50)
    assert not exact
    assert count > 0
