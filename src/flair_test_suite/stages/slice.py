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
from typing import Dict, List

from .base import StageBase             # base class for pipeline stages
from ..lib import PathBuilder          # handles stage directory creation
from ..lib.input_hash import hash_many # hashing inputs for signature
from ..lib.signature import write_marker
from ..lib.reinstate import Reinstate
from ..qc import load_marker            # reuse loader for skip metadata


class SliceStage(StageBase):
    """Stage that slices region-specific subsets of upstream outputs.

    We avoid spawning an external helper script and instead perform slicing
    directly here so that signature / skip / marker logic all remain within
    the existing StageBase framework.
    """
    name = "slice"
    requires = ("align",)
    primary_output_key = "manifest"   # see expected_outputs()

    # ----------------------------- region parsing -----------------------------
    @staticmethod
    def _read_regions(tsv: Path) -> list[tuple[str, int, int]]:
        """Parse a 3-column TSV into (chrom,start,end) tuples.

        Skips blank / commented lines. Emits warnings for malformed lines.
        Raises if no valid regions are parsed.
        """
        regions: list[tuple[str, int, int]] = []
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
                    f"regions TSV line {ln_no}: non-integer coords '{start_s}/{end_s}'",
                    UserWarning,
                )
                continue
            if start_i > end_i:  # tolerate reversed spans
                start_i, end_i = end_i, start_i
            regions.append((chrom, start_i, end_i))
        if not regions:
            raise RuntimeError(f"No valid regions parsed from {tsv}")
        return regions

    # --------------------------- low-level slicers ----------------------------
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
        # Create BAM slice then index
        self._run(
            ["samtools", "view", "-b", str(full_bam), f"{chrom}:{start}-{end}"],
            stdout=open(out_bam, "wb"),
        )
        self._run(["samtools", "index", str(out_bam)])
        # Count reads; emit warning if zero
        try:
            count = subprocess.check_output(
                ["samtools", "view", "-c", str(out_bam)], text=True
            ).strip()
            if count == "0":
                warnings.warn(
                    f"No alignments in {out_bam.name} for region {tag}", UserWarning
                )
        except Exception as e:  # defensive: samtools not installed / etc.
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
        """Slice arbitrary BED-like file using AWK.

        Current predicate performs *full containment* (${col_start}>=start &&
        ${col_end}<=end). To change to *overlap* semantics, adjust the AWK
        condition to `${col_end}>=start && ${col_start}<=end`.
        """
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
        """Slice GTF; retains genes/transcripts fully contained within span.

        This mirrors the legacy slicer logic. Relax by changing the final
        containment checks to any-overlap if desired.
        """
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
            if overlaps
            else (start, end)
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

    # -------------------------- StageBase interface ---------------------------
    def build_cmd(self) -> list[str]:
        """Populate signature pieces.

        We do *not* return a real command because this stage overrides run()
        and performs slicing internally. Returning an empty list keeps
        StageBase code paths happy if they inspect the value.
        """
        cfg = self.cfg
        root = cfg.run.input_root
        data_dir = cfg.run.data_dir

        # Upstream align artifacts
        align_pb: PathBuilder = self.upstreams["align"]
        align_bam = align_pb.stage_dir / f"{self.run_id}_flair.bam"
        align_bed = align_pb.stage_dir / f"{self.run_id}_flair.bed"

        # Choose BED (explicit flag > corrected > align)
        bed_flag = None
        flags_cfg = None
        for st in cfg.run.stages:
            if st.name == "slice":
                flags_cfg = st.flags
                break
        if flags_cfg is None:
            raise RuntimeError("slice stage flags missing in TOML")

        if getattr(flags_cfg, "bed", None):
            bed_flag = (Path(root) / data_dir / flags_cfg.bed).resolve()
        elif "correct" in self.upstreams:  # prefer corrected if available
            corr_pb = self.upstreams["correct"]
            candidate = corr_pb.stage_dir / f"{self.run_id}_all_corrected.bed"
            bed_flag = candidate if candidate.exists() else align_bed
        else:
            bed_flag = align_bed

        # Mandatory flag files
        gtf_path = (Path(root) / data_dir / flags_cfg.gtf).resolve()
        regions_tsv = (Path(root) / data_dir / flags_cfg.regions_tsv).resolve()

        if not gtf_path.exists():
            warnings.warn(f"Slice stage GTF not found: {gtf_path}", UserWarning)
        if not regions_tsv.exists():
            raise RuntimeError(f"regions_tsv file not found: {regions_tsv}")

        # Optional inputs
        opt_keys = [
            "junctions",
            "experiment_5_prime_regions_bed_file",
            "experiment_3_prime_regions_bed_file",
            "reference_5_prime_regions_bed_file",
            "reference_3_prime_regions_bed_file",
        ]
        self._optional_files: dict[str, Path] = {}
        for key in opt_keys:
            val = getattr(flags_cfg, key, None)
            if val:
                p = (Path(root) / data_dir / val).resolve()
                self._optional_files[key] = p
                if not p.exists():
                    warnings.warn(f"Optional slice file missing: {p}", UserWarning)

        # Regions (parse early to fail fast)
        self._regions = self._read_regions(regions_tsv)

        # Inputs to hash for signature reproducibility
        self._hash_inputs = [align_bam, bed_flag, gtf_path, regions_tsv, align_pb.signature]
        if "correct" in self.upstreams:
            self._hash_inputs.append(self.upstreams["correct"].signature)
        for p in self._optional_files.values():
            if p.exists():
                self._hash_inputs.append(p)

        # Store resolved paths for use during run()
        self._align_bam = align_bam
        self._bed_file = bed_flag
        self._gtf_path = gtf_path
        self._regions_tsv = regions_tsv

        # Flags text for signature (purely informational)
        parts = [
            f"--bam={align_bam}",
            f"--bed={bed_flag}",
            f"--gtf={gtf_path}",
            f"--regions-tsv={regions_tsv}",
        ]
        for k, p in self._optional_files.items():
            parts.append(f"--{k.replace('_','-')}={p}")
        self._flags_components = parts
        return []  # placeholder; not executed

    def expected_outputs(self) -> Dict[str, Path]:  # StageBase hook
        return {"manifest": Path("slices_manifest.tsv")}

    def collect_qc(self, pb):  # no QC collector for this stage
        return {}

    # ------------------------------- run() override ---------------------------
    def run(self):  # noqa: C901 (complexity acceptable for orchestrator)
        """Custom run logic: perform slicing internally.

        Mirrors StageBase.run() structure to integrate with skip / marker
        logic while avoiding an external tool invocation.
        """
        import time
        from datetime import datetime, timezone
        import platform

        # Build signature pieces
        self.build_cmd()
        self._input_hashes = hash_many(getattr(self, "_hash_inputs", []))
        self._flags_str = " ".join(getattr(self, "_flags_components", []))
        sig = self.signature

        pb = PathBuilder(self.work_dir, self.run_id, self.name, sig)
        stage_dir = pb.stage_dir
        stage_dir.mkdir(parents=True, exist_ok=True)

        outputs = {k: stage_dir / p for k, p in self.expected_outputs().items()}
        primary = outputs[self.primary_output_key]

        action = Reinstate.decide(
            stage_dir,
            primary,
            needs_qc=False,  # no QC collector registered
            stage_name=self.name,
        )

        if action == "skip":
            print(f"[SKIP] {self.name} complete (sig={sig})")
            # recover previous metadata if present
            marker_f = stage_dir / ".completed.json"
            try:
                meta = load_marker(marker_f)
                pb.metadata = meta.get("qc", {})
            except FileNotFoundError:
                pb.metadata = {}
            return pb

        # (No QC-only path because we have no QC collector)

        # Perform slicing
        start = time.time()
        self._execute_slicing(stage_dir)
        runtime = round(time.time() - start, 2)

        # Write marker
        meta = {
            "stage": self.name,
            "signature": sig,
            "cmd": "(internal slice)",   # no external shell command
            "exit_code": 0,
            "runtime_sec": runtime,
            "started": datetime.fromtimestamp(start, timezone.utc).isoformat(),
            "ended": datetime.now(timezone.utc).isoformat(),
            "host": platform.node(),
            "qc": {},  # placeholder
        }
        write_marker(pb, meta)
        return pb

    # ------------------------- internal slicing driver -----------------------
    def _execute_slicing(self, stage_dir: Path):
        """Iterate user-defined regions and slice all requested inputs."""
        regions_root = stage_dir / "regions"
        regions_root.mkdir(exist_ok=True)

        manifest_rows = []
        for chrom, start, end in self._regions:
            tag = f"{chrom}_{start}_{end}"
            region_dir = regions_root / tag
            region_dir.mkdir(parents=True, exist_ok=True)

            bam_slice = region_dir / f"{tag}.bam"
            bed_slice = region_dir / f"{tag}.bed"
            gtf_slice = region_dir / f"{tag}.gtf"

            # BAM slice
            try:
                self._slice_bam(self._align_bam, chrom, start, end, bam_slice, tag)
            except Exception as e:
                warnings.warn(f"Failed BAM slice {tag}: {e}", UserWarning)

            # BED slice (if bed file exists)
            if self._bed_file and self._bed_file.exists():
                try:
                    self._awk_slice(self._bed_file, chrom, start, end, bed_slice, tag)
                except Exception as e:
                    warnings.warn(f"Failed BED slice {tag}: {e}", UserWarning)

            # GTF slice
            if self._gtf_path.exists():
                try:
                    self._slice_gtf(self._gtf_path, chrom, start, end, gtf_slice, tag)
                except Exception as e:
                    warnings.warn(f"Failed GTF slice {tag}: {e}", UserWarning)

            # Junctions
            junc = self._optional_files.get("junctions")
            if junc and junc.exists():
                self._awk_slice(junc, chrom, start, end, region_dir / f"{tag}.tab", tag)

            # Experiment 5'/3'
            exp5 = self._optional_files.get("experiment_5_prime_regions_bed_file")
            if exp5 and exp5.exists():
                self._awk_slice(exp5, chrom, start, end, region_dir / f"{tag}.exp5.bed", tag)
            exp3 = self._optional_files.get("experiment_3_prime_regions_bed_file")
            if exp3 and exp3.exists():
                self._awk_slice(exp3, chrom, start, end, region_dir / f"{tag}.exp3.bed", tag)

            # Reference 5'/3'
            ref5 = self._optional_files.get("reference_5_prime_regions_bed_file")
            if ref5 and ref5.exists():
                self._awk_slice(ref5, chrom, start, end, region_dir / f"{tag}.ref5.bed", tag)
            ref3 = self._optional_files.get("reference_3_prime_regions_bed_file")
            if ref3 and ref3.exists():
                self._awk_slice(ref3, chrom, start, end, region_dir / f"{tag}.ref3.bed", tag)

            manifest_rows.append([chrom, start, end, str(region_dir)])

        # Write manifest summarising all region directories
        with open(stage_dir / "slices_manifest.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["chrom", "start", "end", "region_dir"])
            w.writerows(manifest_rows)
