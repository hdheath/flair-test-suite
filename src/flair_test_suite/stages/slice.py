# src/flair_test_suite/stages/slice.py
# ------------------------------------
# Slice genomic regions out of upstream FLAIR outputs.
#
# Mandatory flags (in TOML):
#   gtf = "gencode.v48.annotation.gtf"
#   regions_tsv = "regions.tsv"         # 3 columns: chrom  start  end
#
# Optional flags:
#   bed = "override.bed"                # explicit BED to slice (else auto)
#   junctions = "WTC11_all.SJ.out.tab"
#   experiment_5_prime_regions_bed_file = "exp5.bed"
#   experiment_3_prime_regions_bed_file = "exp3.bed"
#   reference_5_prime_regions_bed_file = "ref5.bed"
#   reference_3_prime_regions_bed_file = "ref3.bed"
#
# Upstream usage:
#   - align BAM is always used
#   - BED preference: explicit flag > corrected BED (if correct ran) > align BED
#
# Outputs (under stage_dir):
#   regions/<chrom>_<start>_<end>/*.{bam,bam.bai,bed,gtf,tab,exp5.bed,...}
#   slices_manifest.tsv   (primary output; one line per region)
#
# This stage overrides `run()` to perform internal Python slicing instead of
# calling an external CLI tool.

from __future__ import annotations
import csv
import shlex
import subprocess
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

from .base import StageBase             # base class for pipeline stages
from ..lib import PathBuilder          # handles stage directory creation
from ..lib.input_hash import hash_many # hashing inputs for signature
from ..lib.signature import write_marker
from ..lib.reinstate import Reinstate
from ..qc import load_marker            # reuse loader for skip metadata
from .stage_utils import resolve_path
from ..qc import QC_REGISTRY

class SliceStage(StageBase):
    """Stage that slices region-specific subsets of upstream outputs."""
    name = "slice"
    requires = ("align",)
    primary_output_key = "manifest"

    @staticmethod
    def _read_regions(tsv: Path) -> List[Tuple[str, int, int]]:
        """Parse a 3-column TSV into (chrom,start,end) tuples."""
        regions: List[Tuple[str, int, int]] = []
        for ln_no, raw in enumerate(tsv.read_text().splitlines(), start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                warnings.warn(
                    f"regions TSV line {ln_no}: fewer than 3 columns", UserWarning
                )
                continue
            chrom, start_s, end_s = parts[:3]
            try:
                start_i, end_i = int(start_s), int(end_s)
            except ValueError:
                warnings.warn(
                    f"regions TSV line {ln_no}: non-integer coords '{start_s}/{end_s}'", UserWarning
                )
                continue
            if start_i > end_i:
                start_i, end_i = end_i, start_i
            regions.append((chrom, start_i, end_i))
        if not regions:
            raise RuntimeError(f"No valid regions parsed from {tsv}")
        return regions

    @staticmethod
    def _run(cmd: list[str], **kw):
        """Run a subprocess command, raising on non-zero exit."""
        subprocess.run(cmd, check=True, **kw)

    def _slice_bam(
        self, full_bam: Path, chrom: str, start: int, end: int, out_bam: Path, tag: str
    ):
        """Slice a BAM file using samtools; warn if no reads in slice."""
        if out_bam.exists():
            return
        self._run(
            ["samtools", "view", "-b", str(full_bam), f"{chrom}:{start}-{end}"],
            stdout=open(out_bam, "wb"),
        )
        self._run(["samtools", "index", str(out_bam)])
        try:
            count = subprocess.check_output([
                "samtools", "view", "-c", str(out_bam)
            ], text=True).strip()
            if count == "0":
                warnings.warn(
                    f"No alignments in {out_bam.name} for region {tag}", UserWarning
                )
        except Exception as e:
            warnings.warn(
                f"Failed to count alignments in {out_bam.name}: {e}", UserWarning
            )

    def _awk_slice(
        self,
        in_path: Path,
        chrom: str,
        start: int,
        end: int,
        out_path: Path,
        tag: str,
        col_start=2,
        col_end=3,
    ):
        """Slice arbitrary BED-like file using AWK."""
        if not in_path or not in_path.exists():
            return
        if out_path.exists():
            return
        awk_cmd = (
            f"awk '$1==\"{chrom}\" && ${col_start}>={start} && ${col_end}<={end}' "
            f"{shlex.quote(str(in_path))} > {shlex.quote(str(out_path))}"
        )
        self._run(["bash", "-c", awk_cmd])
        if out_path.exists() and out_path.stat().st_size == 0:
            warnings.warn(f"{out_path.name} is empty for region {tag}", UserWarning)

    def _slice_gtf(
        self, gtf_path: Path, chrom: str, start: int, end: int, out_path: Path, tag: str
    ):
        """Slice GTF; retains genes/transcripts fully contained within span."""
        if out_path.exists():
            return
        import re

        target = chrom if chrom.startswith("chr") else f"chr{chrom}"
        feats, lines, genes = {}, {}, []
        with open(gtf_path) as fh:
            for L in fh:
                if L.startswith("#") or not L.strip():
                    continue
                cols = L.rstrip("\n").split("\t")
                if cols[0] != target:
                    continue
                typ = cols[2]
                st, en = int(cols[3]), int(cols[4])
                if typ == "gene":
                    genes.append((st, en, L))
                m = re.search(r'transcript_id\s+"([^"]+)"', cols[8])
                if not m:
                    continue
                tid = m.group(1)
                feats.setdefault(tid, []).append((st, en))
                lines.setdefault(tid, []).append(L)
        overlaps = [
            iv for ivs in feats.values() for iv in ivs if not (iv[1] < start or iv[0] > end)
        ]
        ns, ne = (
            (min(iv[0] for iv in overlaps), max(iv[1] for iv in overlaps))
            if overlaps else (start, end)
        )
        with open(out_path, "w") as out:
            for g_st, g_en, gL in sorted(genes):
                if g_st >= ns and g_en <= ne:
                    out.write(gL)
            for tid, ivs in feats.items():
                if all(iv[0] >= ns and iv[1] <= ne for iv in ivs):
                    for L in sorted(lines[tid], key=lambda L: int(L.split("\t")[3])):
                        out.write(L)
        if out_path.exists() and out_path.stat().st_size == 0:
            warnings.warn(f"{out_path.name} is empty for region {tag}", UserWarning)

    def build_cmd(self) -> list[str]:
        """Populate signature pieces using shared helpers."""
        cfg = self.cfg
        root = Path(cfg.run.input_root)
        data_dir = Path(cfg.run.data_dir)

        align_pb: PathBuilder = self.upstreams["align"]
        align_bam = align_pb.stage_dir / f"{self.run_id}_flair.bam"
        align_bed = align_pb.stage_dir / f"{self.run_id}_flair.bed"

        flags_cfg = next((st.flags for st in cfg.run.stages if st.name == "slice"), None)
        if flags_cfg is None:
            raise RuntimeError("slice stage flags missing in TOML")

        if getattr(flags_cfg, "bed", None):
            bed_flag = resolve_path(flags_cfg.bed, root=root, data_dir=data_dir)
        elif "correct" in self.upstreams:
            corr_pb = self.upstreams["correct"]
            candidate = corr_pb.stage_dir / f"{self.run_id}_all_corrected.bed"
            bed_flag = candidate if candidate.exists() else align_bed
        else:
            bed_flag = align_bed

        gtf_path = resolve_path(flags_cfg.gtf, root=root, data_dir=data_dir)
        regions_tsv = resolve_path(flags_cfg.regions_tsv, root=root, data_dir=data_dir)
        if not gtf_path.exists():
            warnings.warn(f"Slice stage GTF not found: {gtf_path}", UserWarning)
        if not regions_tsv.exists():
            raise RuntimeError(f"regions_tsv not found: {regions_tsv}")

        opt_keys = [
            "junctions",
            "experiment_5_prime_regions_bed_file",
            "experiment_3_prime_regions_bed_file",
            "reference_5_prime_regions_bed_file",
            "reference_3_prime_regions_bed_file",
        ]
        self._optional_files: Dict[str, Path] = {}
        for key in opt_keys:
            val = getattr(flags_cfg, key, None)
            if val:
                p = resolve_path(val, root=root, data_dir=data_dir)
                self._optional_files[key] = p
                if not p.exists():
                    warnings.warn(f"Optional slice file missing: {p}", UserWarning)

        self._regions = self._read_regions(regions_tsv)
        self._hash_inputs = [align_bam, bed_flag, gtf_path, regions_tsv, align_pb.signature]
        if "correct" in self.upstreams:
            self._hash_inputs.append(self.upstreams["correct"].signature)
        for p in self._optional_files.values():
            if p.exists():
                self._hash_inputs.append(p)

        parts = [
            f"--bam={align_bam}",
            f"--bed={bed_flag}",
            f"--gtf={gtf_path}",
            f"--regions-tsv={regions_tsv}",
        ] + [f"--{k.replace('_','-')}={p}" for k, p in self._optional_files.items()]
        self._flags_components = parts
        # store resolved paths from build_cmd for slicing and QC
        self._align_bam = align_bam
        self._bed_file = bed_flag
        self._gtf_path = gtf_path
        return []

    def expected_outputs(self) -> Dict[str, Path]:
        return {"manifest": Path("slices_manifest.tsv")}  

    def collect_qc(self, pb):
        return {}

    def run(self):
        """Custom run logic: perform slicing internally."""
        import time
        from datetime import datetime, timezone
        import platform

        self.build_cmd()
        self._input_hashes = hash_many(getattr(self, "_hash_inputs", []))
        self._flags_str = " ".join(getattr(self, "_flags_components", []))
        sig = self.signature

        pb = PathBuilder(self.work_dir, self.run_id, self.name, sig)
        stage_dir = pb.stage_dir
        stage_dir.mkdir(parents=True, exist_ok=True)

        manifest = stage_dir / "slices_manifest.tsv"
        start = time.time()
        self._execute_slicing(stage_dir)
        runtime = round(time.time() - start, 2)

        # run region_qc collector if registered
        qc_fn = QC_REGISTRY.get(self.name)
        qc_metrics = {}
        if qc_fn:
            qc_metrics = qc_fn(manifest, stage_dir, runtime_sec=runtime)

        meta = {
            "stage": self.name,
            "signature": sig,
            "cmd": "(internal slice)",
            "exit_code": 0,
            "runtime_sec": runtime,
            "started": datetime.fromtimestamp(start, timezone.utc).isoformat(),
            "ended": datetime.now(timezone.utc).isoformat(),
            "host": platform.node(),
            "qc": qc_metrics,
        }
        write_marker(pb, meta)
        return pb

    def _execute_slicing(self, stage_dir: Path):
        """Iterate user-defined regions and slice all requested inputs."""
        regions_root = stage_dir / "regions"
        regions_root.mkdir(exist_ok=True)

        manifest_rows: List[List[str]] = []
        for chrom, start, end in self._regions:
            tag = f"{chrom}_{start}_{end}"
            region_dir = regions_root / tag
            region_dir.mkdir(parents=True, exist_ok=True)

            bam_slice = region_dir / f"{tag}.bam"
            bed_slice = region_dir / f"{tag}.bed"
            gtf_slice = region_dir / f"{tag}.gtf"

            try:
                self._slice_bam(self._align_bam, chrom, start, end, bam_slice, tag)
            except Exception as e:
                warnings.warn(f"Failed BAM slice {tag}: {e}", UserWarning)

            if self._bed_file and self._bed_file.exists():
                try:
                    self._awk_slice(self._bed_file, chrom, start, end, bed_slice, tag)
                except Exception as e:
                    warnings.warn(f"Failed BED slice {tag}: {e}", UserWarning)

            if self._gtf_path.exists():
                try:
                    self._slice_gtf(self._gtf_path, chrom, start, end, gtf_slice, tag)
                except Exception as e:
                    warnings.warn(f"Failed GTF slice {tag}: {e}", UserWarning)

            for key, fpath in self._optional_files.items():
                out_file = region_dir / f"{tag}.{key.replace('_', '.') }"
                if fpath and fpath.exists():
                    try:
                        self._awk_slice(fpath, chrom, start, end, out_file, tag)
                    except Exception as e:
                        warnings.warn(f"Failed slice for {key} {tag}: {e}", UserWarning)

            manifest_rows.append([chrom, str(start), str(end), str(region_dir)])

        # write manifest
        with open(stage_dir / "slices_manifest.tsv", "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["chrom", "start", "end", "region_dir"] + list(self._optional_files.keys()))
            writer.writerows(manifest_rows)
