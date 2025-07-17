# src/flair_test_suite/qc/align_qc.py
import subprocess
from pathlib import Path
from . import register


@register("align")
def collect(
    bam: Path,
    out_dir: Path,
    runtime_sec: float | None = None,
) -> dict:
    """
    Parse samtools flagstat for total reads, then use the BED line count
    for both mapped_reads and mapped_pct (primary only).
    """
    # 1) total_reads from flagstat
    flag_txt = subprocess.check_output(
        ["samtools", "flagstat", bam], text=True
    )
    total = None
    for line in flag_txt.splitlines():
        if " in total" in line:
            total = int(line.split()[0])
            break
    if total is None:
        raise ValueError("Could not parse total reads from flagstat")

    # 2) BED line count = #primary mapped reads
    bed = bam.with_suffix(".bed")
    if not bed.exists():
        raise FileNotFoundError(f"Expected BED {bed} not found")
    retained = int(subprocess.check_output(["wc", "-l", str(bed)], text=True).split()[0])

    # 3) compute metrics
    mapped_reads   = retained
    mapped_pct     = round(100 * mapped_reads / total, 2) if total else None
    retained_pct   = mapped_pct  # same as mapped_pct in this model

    metrics = {
        "n_total_reads": total,
        "n_mapped_reads": mapped_reads,      # now bed lines only
        "mapped_pct": mapped_pct,            # bed / total
        "n_retained_bed": retained,          # same as mapped_reads
        "retained_pct": retained_pct,        # bed / total
    }
    if runtime_sec is not None:
        metrics["runtime_sec"] = round(runtime_sec, 2)

    # 4) write sidecar TSV
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv = out_dir / "align_qc.tsv"
    with open(tsv, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in metrics.items():
            fh.write(f"{k}\t{v}\n")

    return metrics
