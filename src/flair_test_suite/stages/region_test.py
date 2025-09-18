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
    from ..qc import region_test_qc as _force_import_region_test_qc  # noqa: F401
except Exception:
    pass

__all__ = ["RegionTestStage"]

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


def _parse_tokens(flags_block) -> Dict[str, str | bool]:
    tokens: List[str] = []
    if isinstance(flags_block, list):
        tokens = [str(t).strip() for t in flags_block if str(t).strip()]
    elif isinstance(flags_block, str):
        tokens = [t.strip() for t in flags_block.split(',') if t.strip()]
    else:
        return {}
    out: Dict[str, str | bool] = {}
    import shlex
    for tok in tokens:
        parts = shlex.split(tok)
        if not parts:
            continue
        head = parts[0]
        # --opt=value form
        if '=' in head:
            k, v = head.lstrip('-').split('=', 1)
            out[k.strip()] = v.strip()
            continue
        # --opt value form (or bare presence flag)
        k = head.lstrip('-').strip()
        if len(parts) > 1 and not parts[1].startswith('-'):
            out[k] = parts[1]
        else:
            out[k] = True
    return out


class RegionTestStage(StageBase):
    """
    region_test: materialize per-region BAM/BED/GTF/FASTA + optional inputs.
    Source BED priority: user override via --bed=, else align BED.
    Requires run-level GTF.
    """
    name = "region_test"
    requires = ("align",)
    primary_output_key = "region_details"

    def build_cmds(self) -> List[List[str]]:
        cfg = self.cfg
        data_dir = Path(cfg.run.data_dir)

        # Flags (raw CLI tokens) for this stage only
        stage_cfg = next(st for st in cfg.run.stages if st.name == "region_test")
        flags_map = _parse_tokens(getattr(stage_cfg, "flags", None))
        # Determine overrides from flags only
        bed_override = flags_map.get("bed")
        bam_override = flags_map.get("bam")

        # Upstream align (only required if overrides are not provided)
        align_pb = self.upstreams.get("align")
        if not align_pb and not (bed_override and bam_override):
            raise RuntimeError("region_test requires align upstream or both --bed and --bam")

        if align_pb:
            self._align_bam = align_pb.stage_dir / f"{self.run_id}_flair.bam"
            self._align_bed = align_pb.stage_dir / f"{self.run_id}_flair.bed"
        else:
            # Use user-provided overrides
            self._align_bam = resolve_path(str(bam_override), data_dir=data_dir)
            self._align_bed = resolve_path(str(bed_override), data_dir=data_dir)
            self.logger.info("Using override BAM: %s", self._align_bam)
            self.logger.info("Using override BED: %s", self._align_bed)

        gtf = getattr(cfg.run, "gtf", None)
        if not gtf:
            raise RuntimeError("No GTF specified in run-level inputs for region_test")
        self._gtf_path = resolve_path(gtf, data_dir=data_dir)

        if bed_override:
            self._bed_file = resolve_path(str(bed_override), data_dir=data_dir)
            self.logger.info("Using override BED: %s", self._bed_file)
        elif align_pb and self._align_bed.exists():
            self._bed_file = self._align_bed
            self.logger.info("Using align BED: %s", self._bed_file)
        else:
            raise RuntimeError("No valid BED file found for region_test stage.")

        regions_flag = flags_map.get("regions-tsv") or flags_map.get("regions_tsv")
        if not regions_flag:
            raise RuntimeError("No --regions-tsv specified for region_test stage")
        regions_tsv_path = resolve_path(str(regions_flag), data_dir=data_dir)
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
                self.logger.debug("Optional input: %s -> %s", k, p)

        # Signature inputs
        self._hash_inputs = [
            self._align_bam, self._bed_file, self._gtf_path,
            regions_tsv_path,
            *self._optional.values()
        ]
        if align_pb:
            self._hash_inputs.append(align_pb.signature)
        if self._genome_fa_abs:
            self._hash_inputs.append(Path(self._genome_fa_abs))

        cmds: List[List[str]] = []

        # region_details.tsv (primary) under qc/
        header = "chrom\tstart\tend\tspan_bp"
        body = "\n".join(f"{c}\t{s}\t{e}\t{e - s + 1}" for c, s, e in self._regions)
        cmds.append(["bash", "-lc", "mkdir -p qc"])
        heredoc = f"cat > qc/region_details.tsv << 'EOF'\n{header}\n{body}\nEOF"
        cmds.append(["bash", "-lc", heredoc])

        # Per-region artifacts
        for chrom, start, end in self._regions:
            tag = f"{chrom}_{start}_{end}"
            tmp = f"tmp_sort_{tag}"

            # BAM
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

            # BED (0-based; cols 2-3)
            bed_cmd = (
                "awk -v c='%s' -v s=%d -v e=%d "
                r" -F'\t' '($1==c) && ($2>=s) && ($3<=e){print}' "
                "'%s' > '%s'; test -s '%s' || : > '%s'"
            ) % (chrom, start, end, str(self._bed_file), f"{tag}.bed", f"{tag}.bed", f"{tag}.bed")
            cmds.append(["bash", "-lc", bed_cmd])

            # GTF (1-based; cols 4-5), but restrict to transcripts fully contained in [start,end]
            # Two-step AWK:
            #  1) collect transcript_ids where the transcript feature is fully contained
            #  2) output only features whose transcript_id is in that set and fully contained
            tids_tmp = f"{tag}.tids.tmp"
            awk_tids = (
                "awk -v c='%s' -v s=%d -v e=%d -F'\t' "
                "'($1==c) && ($3==\"transcript\") && ($4>=s) && ($5<=e) {"
                " if (match($9, /transcript_id \"[^\"]+\"/)) {tid=substr($9, RSTART+15, RLENGTH-16); print tid} }' "
                "'%s' | sort -u > '%s'"
            ) % (chrom, start, end, str(self._gtf_path), tids_tmp)
            cmds.append(["bash", "-lc", awk_tids])

            awk_gtf = (
                "awk -v c='%s' -v s=%d -v e=%d -F'\t' "
                "'FNR==NR {a[$1]=1; next} ($1==c) && ($4>=s) && ($5<=e) {"
                " if (match($9, /transcript_id \"[^\"]+\"/)) {tid=substr($9, RSTART+15, RLENGTH-16); if (tid in a) print} }' "
                "'%s' '%s' > '%s'; test -s '%s' || : > '%s'; rm -f '%s'"
            ) % (chrom, start, end, tids_tmp, str(self._gtf_path), f"{tag}.gtf", f"{tag}.gtf", f"{tag}.gtf", tids_tmp)
            cmds.append(["bash", "-lc", awk_gtf])

            # FASTA slice (if genome)
            if self._genome_fa_abs:
                fa_cmd = "samtools faidx '%s' '%s:%d-%d' > '%s' || : > '%s'" % (
                    self._genome_fa_abs, chrom, start, end, f"{tag}.fa", f"{tag}.fa"
                )
                cmds.append(["bash", "-lc", fa_cmd])

            # Optional inputs
            for key, p in self._optional.items():
                if key == "junctions" or str(p).endswith("SJ.out.tab"):
                    dst = f"{tag}.SJ.out.tab"
                    awk = (
                        "awk -v c='%s' -v s=%d -v e=%d "
                        r" -F'\t' '($1==c) && ($2>=s) && ($3<=e){print}' "
                        "'%s' > '%s'; test -s '%s' || : > '%s'"
                    ) % (chrom, start, end, str(p), dst, dst, dst)
                    cmds.append(["bash", "-lc", awk])
                    self.logger.debug("Will slice SJ.out.tab -> %s", dst)
                else:
                    dst = f"{tag}_{Path(p).name}"
                    awk = (
                        "awk -v c='%s' -v s=%d -v e=%d "
                        r" -F'\t' '($1==c) && ($2>=s) && ($3<=e){print}' "
                        "'%s' > '%s'; test -s '%s' || : > '%s'"
                    ) % (chrom, start, end, str(p), dst, dst, dst)
                    cmds.append(["bash", "-lc", awk])
                    self.logger.debug("Will slice BED-like optional -> %s", dst)

        return cmds

    def _region_outputs_present(self, stage_dir: Path) -> bool:
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
            self.logger.error("Command failed (exit %s): %s", proc.returncode, ' '.join(cmd))
            raise RuntimeError(f"{self.name} failed with exit code {proc.returncode}")
        return proc.returncode

    def expected_outputs(self) -> dict[str, Path]:
        return {
            "region_bam": Path("{chrom}_{start}_{end}.bam"),
            "region_bai": Path("{chrom}_{start}_{end}.bam.bai"),
            "region_bed": Path("{chrom}_{start}_{end}.bed"),
            "region_gtf": Path("{chrom}_{start}_{end}.gtf"),
            "region_fa": Path("{chrom}_{start}_{end}.fa"),
            "region_details": Path("qc/region_details.tsv"),
        }
