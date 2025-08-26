#!/usr/bin/env python3
"""SQANTI QC collector for collapse/transcriptome stages."""

from __future__ import annotations

import csv
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from . import register, write_metrics
from ..stages import stage_utils

logger = logging.getLogger(__name__)

# Optional plotting dependency
try:  # pragma: no cover - heavy deps optional
    from ..plotting.sqanti_plot import plot_summary  # type: ignore
except Exception:  # pragma: no cover - handle missing deps
    plot_summary = None  # type: ignore

# ---------------------------------------------------------------------------
# SQANTI classification parsing helpers
ISOCLASSES = [
    "FSM",
    "ISM",
    "NIC",
    "NNC",
    "Genic_Genomic",
    "Genic_Intron",
    "Antisense",
    "Fusion",
    "Intergenic",
]
_raw_keys = {
    "full-splice_match": "FSM",
    "incomplete-splice_match": "ISM",
    "novel_in_catalog": "NIC",
    "novel-not-in-catalog": "NNC",
    "novel_not_in_catalog": "NNC",
    "genic": "Genic_Genomic",
    "genic_intron": "Genic_Intron",
    "genic_genomic": "Genic_Genomic",
    "antisense": "Antisense",
    "fusion": "Fusion",
    "intergenic": "Intergenic",
}
NORMALIZED_MAP = {k.strip().lower().replace("-", "_"): v for k, v in _raw_keys.items()}


def _conda_env_exists(env: str) -> bool:
    """Return True if the named Conda environment can be activated."""
    try:
        subprocess.run(
            ["conda", "run", "-n", env, "python", "-V"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        return True
    except Exception:
        return False


def _run_sqanti(
    gtf: Path,
    ref_gtf: Path,
    ref_fa: Path,
    prime5: Optional[Path],
    prime3: Optional[Path],
    sj: Optional[Path],
    prefix: str,
    env: str,
    cpus: int,
    out_dir: Path,
    log_dir: Path,
) -> None:
    # Prefer invoking the package via `python -m` inside the conda env to avoid
    # relying on an installed entrypoint script. If that fails with command not
    # found, fall back to running the entrypoint through a shell.
    cmd = [
        "conda",
        "run",
        "-n",
        env,
        "sqanti3_qc.py",
        str(gtf),
        str(ref_gtf),
        str(ref_fa),
    ]
    if prime5:
        cmd.extend(["--CAGE_peak", str(prime5)])
    if prime3:
        cmd.extend(["--polyA_peak", str(prime3)])
    if sj:
        cmd.extend(["-c", str(sj)])
    # sqanti3_qc uses -t for number of threads (CPUs)
    cmd.extend(["-o", prefix, "-d", str(out_dir), "--skipORF", "-t", str(cpus)])

    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "sqanti.log"
    # run and capture output; if returncode indicates command not found (127),
    # try a shell fallback. Surface helpful log excerpt on failure.
    try:
        with open(log_path, "w") as lf:
            subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, check=True)
    except subprocess.CalledProcessError as e:
        # If the python -m invocation failed with code 127, try a shell entrypoint
        if e.returncode == 127:
            logger.warning(f"[SQANTI] python -m invocation failed (rc=127); attempting shell fallback for env '{env}'")
            shell_cmd = (
                "conda run -n " + env + " bash -lc 'sqanti3_qc.py "
                + "\"" + str(gtf) + "\" "
                + "\"" + str(ref_gtf) + "\" "
                + "\"" + str(ref_fa) + "\""
            )
            if sj:
                shell_cmd += " -c \"%s\"" % str(sj)
            shell_cmd += " -o \"%s\" -d \"%s\" --skipORF -t %s'" % (prefix, str(gtf.parent), str(cpus))
            try:
                with open(log_path, "a") as lf:
                    subprocess.run(shell_cmd, stdout=lf, stderr=subprocess.STDOUT, shell=True, check=True)
                return
            except subprocess.CalledProcessError:
                logger.warning(f"[SQANTI] shell fallback also failed for env '{env}'")
        # For any failure, include a short excerpt of the log to help debugging
        try:
            with open(log_path, "r") as lf:
                lines = lf.readlines()[-200:]
                excerpt = "".join(lines[-50:]) if lines else "(no log content)"
        except Exception:
            excerpt = "(could not read sqanti.log)"
        logger.warning(f"[SQANTI] sqanti3_qc failed for {prefix}: {e}; log excerpt:\n{excerpt}")
        # Re-raise so callers can still observe the failure if desired
        raise


def _normalize(lbl: str) -> str:
    return lbl.strip().lower().replace("-", "_")


def _safe_float(v: str) -> float:
    try:
        return float(v)
    except Exception:
        return float("nan")


def _summarize(class_files: List[Path], out_tsv: Path) -> None:
    """Summarize SQANTI classification files into ``out_tsv``."""
    rows: List[List[object]] = []
    for path in class_files:
        if not path.exists():
            logger.warning(f"[SQANTI] missing classification file: {path}")
            continue
        tot = cages = polya = ref5 = ref3 = srtm = sntm = 0
        counts: Dict[str, int] = {c: 0 for c in ISOCLASSES}
        warned: set[str] = set()
        with path.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if cols[0] == "isoform" or len(cols) < 44:
                    continue
                tot += 1
                raw = cols[5]
                iso = NORMALIZED_MAP.get(_normalize(raw))
                if iso in counts:
                    counts[iso] += 1
                else:
                    if raw not in warned:
                        logger.warning(f"[SQANTI] Unknown class '{raw}' in {path}")
                        warned.add(raw)
                dC = _safe_float(cols[39]); iC = cols[40].lower() == "true"
                dP = _safe_float(cols[41]); iP = cols[42].lower() == "true"
                hasC = iC or (not (dC != dC) and abs(dC) <= 50)
                hasP = iP or (not (dP != dP) and abs(dP) <= 50)
                dTSS = _safe_float(cols[10]); dTTS = _safe_float(cols[11])
                gTSS = _safe_float(cols[12]); gTTS = _safe_float(cols[13])
                a5 = ((not (dTSS != dTSS) and abs(dTSS) <= 50) or (not (gTSS != gTSS) and abs(gTSS) <= 50))
                a3 = ((not (dTTS != dTTS) and abs(dTTS) <= 50) or (not (gTTS != gTTS) and abs(gTTS) <= 50))
                if a5:
                    ref5 += 1
                if a3:
                    ref3 += 1
                if iso in ("FSM", "ISM"):
                    if (hasC or a5) and (hasP or a3):
                        srtm += 1
                else:
                    if (hasC or a5) and (hasP or a3):
                        sntm += 1
                cages += hasC
                polya += hasP
        sample = path.name.replace("_classification.txt", "")
        rows.append([
            sample,
            round(srtm / tot, 3) if tot else 0,
            round(sntm / tot, 3) if tot else 0,
            round(cages / tot, 3) if tot else 0,
            round(polya / tot, 3) if tot else 0,
            round(ref5 / tot, 3) if tot else 0,
            round(ref3 / tot, 3) if tot else 0,
            *[counts[c] for c in ISOCLASSES],
            tot,
        ])
    header = [
        "sample",
        "fracSRTM",
        "fracSNTM",
        "5'endsupport",
        "3'endsupport",
        "5'refsupport",
        "3'refsupport",
        *ISOCLASSES,
        "total_isoforms",
    ]
    with out_tsv.open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    logger.info(f"[SQANTI] Wrote summary {out_tsv}")


# ---------------------------------------------------------------------------
# Collector entry point
@register("SQANTI")
@register("sqanti")
@register("collapse")
@register("transcriptome")
def collect(stage_dir: Path, cfg, upstreams=None, *, out_dir: Optional[Path] = None) -> None:
    """Run SQANTI QC on all ``*.isoforms.gtf`` under ``stage_dir``.

    Parameters
    ----------
    stage_dir : Path
        Directory containing isoform GTFs.
    cfg : object
        Config object.
    out_dir : Path | None
        Directory to write SQANTI outputs. Defaults to ``stage_dir``.
    """
    stage_name = stage_dir.parent.name
    qc_block = getattr(cfg, "qc", {}).get(stage_name, {}).get("SQANTI", {})
    env = qc_block.get("conda_env") if isinstance(qc_block, dict) else None
    if not env:
        logger.info("[SQANTI] conda_env not configured; skipping")
        return
    if not _conda_env_exists(env):
        logger.warning(f"[SQANTI] conda env '{env}' not found; skipping")
        return
    cpus = int(qc_block.get("cpus", 1)) if isinstance(qc_block, dict) else 1

    data_dir = Path(cfg.run.data_dir)
    ref_gtf = stage_utils.resolve_path(cfg.run.gtf, data_dir=data_dir)
    ref_fa = stage_utils.resolve_path(cfg.run.genome_fa, data_dir=data_dir)
    prime5 = getattr(cfg.run, "experiment_5_prime_regions_bed_file", None)
    prime3 = getattr(cfg.run, "experiment_3_prime_regions_bed_file", None)
    sj = getattr(cfg.run, "junctions", None)
    prime5 = stage_utils.resolve_path(prime5, data_dir=data_dir) if prime5 else None
    prime3 = stage_utils.resolve_path(prime3, data_dir=data_dir) if prime3 else None
    sj = stage_utils.resolve_path(sj, data_dir=data_dir) if sj else None

    gtfs = sorted(stage_dir.glob("*.isoforms.gtf"))
    if not gtfs:
        logger.warning("[SQANTI] no isoforms.gtf found under %s", stage_dir)
        return

    out_dir = out_dir or stage_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    logs_base = out_dir / "logs"
    class_files: List[Path] = []
    for gtf in gtfs:
        prefix = gtf.name.replace(".isoforms.gtf", "")
        class_txt = out_dir / f"{prefix}_classification.txt"
        class_files.append(class_txt)
        if class_txt.exists():
            logger.info(f"[SQANTI] skipping SQANTI for {gtf.name}; classification present")
            continue
        try:
            _run_sqanti(
                gtf,
                ref_gtf,
                ref_fa,
                prime5,
                prime3,
                sj,
                prefix,
                env,
                cpus,
                out_dir,
                logs_base / prefix,
            )
        except Exception as e:
            logger.warning(f"[SQANTI] sqanti3_qc failed for {gtf.name}: {e}")

    existing = [p for p in class_files if p.exists()]
    if not existing:
        logger.warning("[SQANTI] no classification files found; skipping summary")
        return
    out_tsv = out_dir / "sqanti_results.tsv"
    _summarize(existing, out_tsv)
    if plot_summary:
        try:  # pragma: no cover - plotting is heavy
            plot_summary(str(out_tsv), out_dir)
        except Exception as e:  # pragma: no cover - best effort
            logger.warning(f"[SQANTI] plot failed: {e}")
    write_metrics(out_dir, "SQANTI", {"tsv": str(out_tsv)})
    logger.info("[SQANTI] QC complete")
