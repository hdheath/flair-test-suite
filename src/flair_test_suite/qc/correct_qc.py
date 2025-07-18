from __future__ import annotations
import subprocess
from pathlib import Path
from . import register


@register("correct")
def collect(bed: Path, out_dir: Path,
            n_input_reads: int,
            runtime_sec: float | None = None) -> dict:

    retained = int(subprocess.check_output(["wc", "-l", str(bed)], text=True).split()[0])
    retained_pct = (round(100 * retained / n_input_reads, 2)
                if n_input_reads else None)


    metrics = {
        "n_input_reads":     n_input_reads,
        "n_corrected_reads": retained,
        "retained_pct":      retained_pct,
    }
    if runtime_sec is not None:
        metrics["runtime_sec"] = round(runtime_sec, 2)

    (out_dir / "correct_qc.tsv").write_text(
        "metric\tvalue\n" + "\n".join(f"{k}\t{v}" for k, v in metrics.items())
    )
    return metrics
