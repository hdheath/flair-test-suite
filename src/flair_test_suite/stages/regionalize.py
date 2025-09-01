from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Tuple

from .base import StageBase, StageAction
from .stage_utils import resolve_path
from ..lib.paths import PathBuilder
from ..lib.input_hash import hash_many
import subprocess


logger = logging.getLogger(__name__)

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
            logger.warning("regions TSV line %d: fewer than 3 columns", ln)
            continue
        try:
            chrom = parts[0]; s = int(parts[1]); e = int(parts[2])
        except ValueError:
            logger.warning("regions TSV line %d: non-integer coords", ln)
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

        gtf = getattr(cfg.run, "gtf", None) or flags.get("gtf")
        if not gtf:
            raise RuntimeError("No GTF specified in regionalize inputs")
        self._gtf_path = resolve_path(gtf, data_dir=data_dir)

        override_bed = flags.get("bed")
        if override_bed:
            self._bed_file = resolve_path(override_bed, data_dir=data_dir)
            self.logger.info("Using override BED: %s", self._bed_file)
        elif self._align_bed.exists():
            self._bed_file = self._align_bed
            self.logger.info("Using align BED: %s", self._bed_file)
        else:
            raise RuntimeError("No valid BED file found for regionalize stage.")

        regions_tsv = getattr(cfg.run, "regions_tsv", None) or flags.get("regions_tsv")
        if not regions_tsv:
            raise RuntimeError("No regions_tsv specified in regionalize inputs")
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

        # Optional inputs (slice per region) sourced solely from run-level config
        opt_keys = [
            "junctions",  # STAR SJ.out.tab
            "experiment_5_prime_regions_bed_file",
            "experiment_3_prime_regions_bed_file",
            "reference_5_prime_regions_bed_file",
            "reference_3_prime_regions_bed_file",
        ]
        self._optional: Dict[str, Path] = {}
        for k in opt_keys:
            v = getattr(cfg.run, k, None)
            if v:
                p = resolve_path(v, data_dir=data_dir)
                self._optional[k] = p
                self.logger.info("Optional input: %s -> %s", k, p)

        # Signature inputs
        self._hash_inputs = [
            self._align_bam, self._bed_file, self._gtf_path,
            regions_tsv_path, align_pb.signature, *self._optional.values()
        ]
        if self._genome_fa_abs:
            self._hash_inputs.append(Path(self._genome_fa_abs))

        cmds: List[List[str]] = []

        # region_details.tsv (primary) under qc/ (avoid repeating 'regionalize')
        header = "chrom\tstart\tend\tspan_bp"
        body = "\n".join(f"{c}\t{s}\t{e}\t{e - s + 1}" for c, s, e in self._regions)
        cmds.append(["bash", "-lc", "mkdir -p qc"])
        heredoc = f"cat > qc/region_details.tsv << 'EOF'\n{header}\n{body}\nEOF"
        cmds.append(["bash", "-lc", heredoc])

        # Per-region artifacts
        for chrom, start, end in self._regions:
            tag = f"{chrom}_{start}_{end}"
            tmp = f"tmp_sort_{tag}"

            # BAM (fully contained reads only):
            # Use samtools view -h to get SAM, filter with awk to keep only
            # alignments whose reference span is fully within [start,end], then
            # convert back to BAM and sort.
            bam_pipe = (
                "samtools view -h '%(bam)s' '%(chrom)s:%(start)d-%(end)d' | "
                "awk -v s=%(start)d -v e=%(end)d '"
                "BEGIN{OFS=\"\t\"} "
                "/^@/ {print; next} "
                "{pos=$4; cig=$6; if (pos==\"\" || cig==\"*\") next; "
                "len=0; c=cig; while (match(c, /[0-9]+[MIDNSHP=X]/)) {n=substr(c, RSTART, RLENGTH-1); op=substr(c, RSTART+RLENGTH-1, 1); if (op ~ /[MDN=X]/) len += n; c=substr(c, RSTART+RLENGTH);} "
                "end=pos+len-1; if (pos>=s && end<=e) print}' | "
                "samtools view -Sb - | samtools sort -o '%(tag)s.bam' -T '%(tmp)s'"
            ) % {"bam": str(self._align_bam), "chrom": chrom, "start": start, "end": end, "tag": tag, "tmp": tmp}
            cmds.append(["bash", "-lc", bam_pipe])
            cmds.append(["bash", "-lc", f"samtools index '{tag}.bam' || :"])

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
                    self.logger.info("Will slice SJ.out.tab -> %s", dst)
                else:
                    # Treat as BED-like (0-based, cols 2-3)
                    dst = f"{tag}_{Path(p).name}"
                    awk = (
                        "awk -v c='%s' -v s=%d -v e=%d "
                        r" -F'\t' '($1==c) && ($2>=s) && ($3<=e){print}' "
                        "'%s' > '%s'; test -s '%s' || : > '%s'"
                    ) % (chrom, start, end, str(p), dst, dst, dst)
                    cmds.append(["bash", "-lc", awk])
                    self.logger.info("Will slice BED-like optional -> %s", dst)

        return cmds

    def _region_outputs_present(self, stage_dir: Path) -> bool:
        """Return True if all per-region artifacts exist and are non-empty.

        Checks that for every region in qc/regionalize/region_details.tsv, the
        corresponding BAM and GTF files are present (and non-empty). BED is
        optional for this completeness test.
        """
        details = stage_dir / "qc" / "region_details.tsv"
        if not details.exists():
            return False
        try:
            regions = []
            with open(details) as fh:
                next(fh, None)
                for raw in fh:
                    line = raw.strip()
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) < 3:
                        continue
                    chrom, s, e = parts[0], parts[1], parts[2]
                    regions.append(f"{chrom}_{s}_{e}")
            if not regions:
                return False
            for tag in regions:
                bam = stage_dir / f"{tag}.bam"
                gtf = stage_dir / f"{tag}.gtf"
                if not (bam.exists() and bam.stat().st_size > 0 and gtf.exists() and gtf.stat().st_size > 0):
                    return False
            return True
        except Exception:
            return False

    def _decide_action(self, stage_dir: Path, primary: Path, needs_qc: bool) -> StageAction:
        """Override to force full run when region outputs are missing.

        Reinstate may suggest QC-only when the primary exists but QC is missing.
        For regionalize, we also require that per-region BAM/GTF files exist;
        if they don't, force a full run to materialize them.
        """
        decision = super()._decide_action(stage_dir, primary, needs_qc)
        if decision == StageAction.QC_ONLY and not self._region_outputs_present(stage_dir):
            self.logger.info("Per-region artifacts missing; forcing full run")
            return StageAction.RUN
        return decision

    def run_tool(
        self,
        cmd: list[str],
        log_path: Path | None = None,
        cwd: Path | None = None,
    ) -> int:
        """Run commands exactly as given (no conda), logging to tool.log."""

        log_path = log_path or Path("tool.log")
        log_path.parent.mkdir(parents=True, exist_ok=True)

        self.logger.info("Running: %s", ' '.join(cmd))
        with open(log_path, "a") as logf:
            proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=cwd)

        if proc.returncode != 0:
            self.logger.error(
                "Command failed (exit %s): %s", proc.returncode, ' '.join(cmd)
            )
            raise RuntimeError(f"{self.name} failed with exit code {proc.returncode}")
        return proc.returncode


    # for legacy callers
    def expected_outputs(self) -> dict[str, Path]:
        return {
            "region_bam": Path("{chrom}_{start}_{end}.bam"),
            "region_bai": Path("{chrom}_{start}_{end}.bam.bai"),
            "region_bed": Path("{chrom}_{start}_{end}.bed"),
            "region_gtf": Path("{chrom}_{start}_{end}.gtf"),
            "region_fa": Path("{chrom}_{start}_{end}.fa"),
            # Primary now lives under qc/
            "region_details": Path("qc/region_details.tsv"),
        }
