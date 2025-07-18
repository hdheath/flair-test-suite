# src/flair_test_suite/qc/align_qc.py
from __future__ import annotations
import subprocess
from pathlib import Path
from . import register


@register("align")
def collect(
    bam: Path,
    out_dir: Path,
    n_input_reads: int,        # StageBase passes this
    runtime_sec: float | None = None,
) -> dict:
    """
    QC for flair align:
      • n_total_reads  – sequences fed to FLAIR (counted pre‑run)
      • n_retained_bed – primary alignments in BED
      • mapped_pct     – retained / total ×100
    """
    bed = bam.with_suffix(".bed")
    if not bed.exists():
        raise FileNotFoundError(bed)

    retained = int(
        subprocess.check_output(["wc", "-l", str(bed)], text=True).split()[0]
    )

    mapped_pct   = round(100 * retained / n_input_reads, 2)
    retained_pct = mapped_pct        # identical definition

    metrics = {
        "n_total_reads":  n_input_reads,
        "n_retained_bed": retained,
        "mapped_pct":     mapped_pct,
        "retained_pct":   retained_pct,
    }
    if runtime_sec is not None:
        metrics["runtime_sec"] = round(runtime_sec, 2)

    # -------- write side‑car TSV ----------
    tsv = out_dir / "align_qc.tsv"
    with tsv.open("w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in metrics.items():
            fh.write(f"{k}\t{v}\n")

    return metrics
