from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Tuple

from .base import StageBase
from .stage_utils import resolve_path
from ..lib.paths import PathBuilder
from ..lib.input_hash import hash_many
import subprocess

# Force-load QC so it's registered
try:
    from ..qc import regionalize_qc as _force_import_regionalize_qc  # noqa: F401
except Exception:
    pass

__all__ = ["RegionalizeStage"]

Region = Tuple[str, int, int]


def _read_regions(tsv: Path) -> List[Region]:
    regions: List[Region] = []
    for ln, raw in enumerate(tsv.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            logging.warning(f"[regionalize] regions TSV line {ln}: fewer than 3 columns")
            continue
        try:
            chrom = parts[0]; s = int(parts[1]); e = int(parts[2])
        except ValueError:
            logging.warning(f"[regionalize] regions TSV line {ln}: non-integer coords")
            continue
        if s > e:
            s, e = e, s
        regions.append((chrom, s, e))
    if not regions:
        raise RuntimeError(f"No valid regions parsed from {tsv}")
    return regions


class RegionalizeStage(StageBase):
    """
    Regionalize: materialize per-region BAM/BED/GTF/FASTA + optional inputs.
    Source BED priority: user override via flags.bed, else align BED.
    """
    name = "regionalize"
    requires = ("align",)
    primary_output_key = "region_details"

    def build_cmds(self) -> List[List[str]]:
        cfg = self.cfg
        data_dir = Path(cfg.run.data_dir)

        # Upstream align
        align_pb = self.upstreams.get("align")
        if not align_pb:
            raise RuntimeError("regionalize requires align upstream")
        self._align_bam = align_pb.stage_dir / f"{self.run_id}_flair.bam"
        self._align_bed = align_pb.stage_dir / f"{self.run_id}_flair.bed"

        # Flags
        flags = next(st.flags for st in cfg.run.stages if st.name == "regionalize")

        gtf = flags.get("gtf")
        if not gtf:
            raise RuntimeError("No GTF specified in regionalize stage flags.")
        self._gtf_path = resolve_path(gtf, data_dir=data_dir)

        override_bed = flags.get("bed")
        if override_bed:
            self._bed_file = resolve_path(override_bed, data_dir=data_dir)
            logging.info(f"[regionalize] Using override BED: {self._bed_file}")
        elif self._align_bed.exists():
            self._bed_file = self._align_bed
            logging.info(f"[regionalize] Using align BED: {self._bed_file}")
        else:
            raise RuntimeError("No valid BED file found for regionalize stage.")

        regions_tsv = flags.get("regions_tsv")
        if not regions_tsv:
            raise RuntimeError("No regions_tsv specified in regionalize stage flags.")
        regions_tsv_path = resolve_path(regions_tsv, data_dir=data_dir)
        self._regions: List[Region] = _read_regions(regions_tsv_path)

        # Optional genome for FASTA slice
        self._genome_fa_abs = None
        genome_fa = getattr(cfg.run, "genome_fa", None)
        if genome_fa:
            try:
                self._genome_fa_abs = str(resolve_path(genome_fa, data_dir=data_dir))
            except Exception:
                self._genome_fa_abs = None

        # Optional inputs (slice per region)
        opt_keys = [
            "junctions",  # STAR SJ.out.tab
            "experiment_5_prime_regions_bed_file",
            "experiment_3_prime_regions_bed_file",
            "reference_5_prime_regions_bed_file",
            "reference_3_prime_regions_bed_file",
        ]
        self._optional: Dict[str, Path] = {}
        for k in opt_keys:
            v = flags.get(k)
            if v:
                p = resolve_path(v, data_dir=data_dir)
                self._optional[k] = p
                logging.info(f"[regionalize] Optional input: {k} -> {p}")

        # Signature inputs
        self._hash_inputs = [
            self._align_bam, self._bed_file, self._gtf_path,
            regions_tsv_path, align_pb.signature, *self._optional.values()
        ]
        if self._genome_fa_abs:
            self._hash_inputs.append(Path(self._genome_fa_abs))

        cmds: List[List[str]] = []

        # region_details.tsv (primary)
        header = "chrom\tstart\tend\tspan_bp"
        body = "\n".join(f"{c}\t{s}\t{e}\t{e - s + 1}" for c, s, e in self._regions)
        heredoc = f"cat > region_details.tsv << 'EOF'\n{header}\n{body}\nEOF"
        cmds.append(["bash", "-lc", heredoc])

        # Per-region artifacts
        for chrom, start, end in self._regions:
            tag = f"{chrom}_{start}_{end}"
            tmp = f"tmp_sort_{tag}"

            # BAM
            pipe = (
                f"samtools view -b '{self._align_bam}' '{chrom}:{start}-{end}'"
                f" | samtools sort -o '{tag}.bam' -T '{tmp}'"
            )
            cmds.append(["bash", "-lc", pipe])
            cmds.append(["bash", "-lc", f"samtools index '{tag}.bam'"])

            # BED (assume 0-based; cols 2-3)
            bed_cmd = (
                "awk -v c='%s' -v s=%d -v e=%d "
                r" -F'\t' '($1==c) && ($2>=s) && ($3<=e){print}' "
                "'%s' > '%s'; test -s '%s' || : > '%s'"
            ) % (chrom, start, end, str(self._bed_file), f"{tag}.bed", f"{tag}.bed", f"{tag}.bed")
            cmds.append(["bash", "-lc", bed_cmd])

            # GTF (1-based; cols 4-5)
            gtf_cmd = (
                "awk -v c='%s' -v s=%d -v e=%d "
                r" -F'\t' '($1==c) && ($4>=s) && ($5<=e){print}' "
                "'%s' > '%s'; test -s '%s' || : > '%s'"
            ) % (chrom, start, end, str(self._gtf_path), f"{tag}.gtf", f"{tag}.gtf", f"{tag}.gtf")
            cmds.append(["bash", "-lc", gtf_cmd])

            # FASTA slice (if genome)
            if self._genome_fa_abs:
                fa_cmd = "samtools faidx '%s' '%s:%d-%d' > '%s' || : > '%s'" % (
                    self._genome_fa_abs, chrom, start, end, f"{tag}.fa", f"{tag}.fa"
                )
                cmds.append(["bash", "-lc", fa_cmd])

            # ----- OPTIONAL INPUTS -----
            for key, p in self._optional.items():
                # STAR junctions: SJ.out.tab (1-based positions in cols 2-3)
                if key == "junctions" or str(p).endswith("SJ.out.tab"):
                    dst = f"{tag}.SJ.out.tab"
                    awk = (
                        "awk -v c='%s' -v s=%d -v e=%d "
                        r" -F'\t' '($1==c) && ($2>=s) && ($3<=e){print}' "
                        "'%s' > '%s'; test -s '%s' || : > '%s'"
                    ) % (chrom, start, end, str(p), dst, dst, dst)
                    cmds.append(["bash", "-lc", awk])
                    logging.info(f"[regionalize] will slice SJ.out.tab -> {dst}")
                else:
                    # Treat as BED-like (0-based, cols 2-3)
                    dst = f"{tag}_{Path(p).name}"
                    awk = (
                        "awk -v c='%s' -v s=%d -v e=%d "
                        r" -F'\t' '($1==c) && ($2>=s) && ($3<=e){print}' "
                        "'%s' > '%s'; test -s '%s' || : > '%s'"
                    ) % (chrom, start, end, str(p), dst, dst, dst)
                    cmds.append(["bash", "-lc", awk])
                    logging.info(f"[regionalize] will slice BED-like optional -> {dst}")

        return cmds

    def run_tool(self, cmd: list[str], log_path: Path | None = None, cwd: Path | None = None) -> int:
        """Run commands exactly as given (no conda), logging to tool.log."""

        log_path = log_path or Path("tool.log")
        log_path.parent.mkdir(parents=True, exist_ok=True)

        logging.info(f"[{self.name}] Running: {' '.join(cmd)}")
        with open(log_path, "a") as logf:
            proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=cwd)

        if proc.returncode != 0:
            logging.error(f"[{self.name}] Command failed (exit {proc.returncode}): {' '.join(cmd)}")
            raise RuntimeError(f"{self.name} failed with exit code {proc.returncode}")
        return proc.returncode


        logging.info(f"[{self.name}] Running (no conda): {' '.join(cmd)}")
        with open(log_path, "a") as logf:
            proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=cwd)

        if proc.returncode != 0:
            logging.error(f"[{self.name}] Command failed (exit {proc.returncode}): {' '.join(cmd)}")
            raise RuntimeError(f"{self.name} failed with exit code {proc.returncode}")
        return proc.returncode

    # for legacy callers
    def build_cmd(self) -> List[str]:
        cmds = self.build_cmds()
        return cmds[0] if cmds else []

    def expected_outputs(self) -> dict[str, Path]:
        return {
            "region_bam": Path("{chrom}_{start}_{end}.bam"),
            "region_bai": Path("{chrom}_{start}_{end}.bam.bai"),
            "region_bed": Path("{chrom}_{start}_{end}.bed"),
            "region_gtf": Path("{chrom}_{start}_{end}.gtf"),
            "region_fa": Path("{chrom}_{start}_{end}.fa"),
            "region_details": Path("region_details.tsv"),
        }
