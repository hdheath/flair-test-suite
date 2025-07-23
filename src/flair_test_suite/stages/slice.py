# src/flair_test_suite/stages/slice.py
# -------------------------------------------------
# Combined-output slicer (v3) – robust region containment
#   • Emits one combined BAM, BED, GTF, and any optional extras
#     for all specified regions.
#   • Writes slice_region_details.tsv (chrom, start, end, span_bp).
#   • Optional .bed and .tab files parsed with correct coords.

from __future__ import annotations
from pathlib import Path
import csv, subprocess, warnings, shlex
from collections import defaultdict
from typing import Dict, List, Tuple

from .base import StageBase
from ..lib import PathBuilder
from ..lib.input_hash import hash_many
from ..lib.signature import write_marker
from ..lib.reinstate import Reinstate
from .stage_utils import resolve_path
from ..qc import QC_REGISTRY
from .stage_utils import resolve_path, filter_file_by_regions

__all__ = ["SliceStage"]

Region = Tuple[str, int, int]  # chrom, start, end (1-based inclusive)


def _read_regions(tsv: Path) -> List[Region]:
    """Parse TSV → List[(chrom,start,end)]."""
    regions: List[Region] = []
    for ln, raw in enumerate(tsv.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            warnings.warn(f"regions TSV line {ln}: fewer than 3 columns", UserWarning)
            continue
        try:
            chrom = parts[0]
            s = int(parts[1]); e = int(parts[2])
        except ValueError:
            warnings.warn(f"regions TSV line {ln}: non-integer coords", UserWarning)
            continue
        if s > e:  # swap
            s, e = e, s
        regions.append((chrom, s, e))
    if not regions:
        raise RuntimeError(f"No valid regions parsed from {tsv}")
    return regions


class IntervalLookup:
    """Efficient full-containment lookup by chromosome."""
    def __init__(self, regions: List[Region]):
        idx: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
        for c, s, e in regions:
            idx[c].append((s, e))
        for c, lst in idx.items():
            lst.sort()
        self.idx = idx

    def contains(self, chrom: str, start: int, end: int) -> bool:
        """Return True if [start,end] fully within any region on chrom."""
        if chrom not in self.idx:
            return False
        for s, e in self.idx[chrom]:
            if start >= s and end <= e:
                return True
            if start < s:
                break
        return False


class SliceStage(StageBase):
    """Stage: combine region slices into unified outputs."""
    name = "slice"
    requires = ("align",)
    primary_output_key = "combined_bam"

    def build_cmd(self):
        cfg = self.cfg
        root = Path(cfg.run.input_root)
        data_dir = Path(cfg.run.data_dir)

        # upstream align outputs
        align_pb = self.upstreams["align"]
        self._align_bam = align_pb.stage_dir / f"{self.run_id}_flair.bam"
        self._align_bed = align_pb.stage_dir / f"{self.run_id}_flair.bed"

        flags = next(st.flags for st in cfg.run.stages if st.name == "slice")

        # select BED file
        if getattr(flags, "bed", None):
            self._bed_file = resolve_path(flags.bed, root=root, data_dir=data_dir)
        elif "correct" in self.upstreams:
            cand = self.upstreams["correct"].stage_dir / f"{self.run_id}_all_corrected.bed"
            self._bed_file = cand if cand.exists() else self._align_bed
        else:
            self._bed_file = self._align_bed

        # mandatory GTF + regions tsv
        self._gtf_path = resolve_path(flags.gtf, root=root, data_dir=data_dir)
        regions_tsv    = resolve_path(flags.regions_tsv, root=root, data_dir=data_dir)
        if not regions_tsv.exists():
            raise RuntimeError(f"regions_tsv not found: {regions_tsv}")

        # parse regions
        self._regions   = _read_regions(regions_tsv)
        self._intervals = [f"{c}:{s}-{e}" for c, s, e in self._regions]
        self._lookup    = IntervalLookup(self._regions)

        # optional extras
        opt_keys = [
            "junctions",
            "experiment_5_prime_regions_bed_file",
            "experiment_3_prime_regions_bed_file",
            "reference_5_prime_regions_bed_file",
            "reference_3_prime_regions_bed_file",
        ]
        self._optional: Dict[str, Path] = {}
        for k in opt_keys:
            v = getattr(flags, k, None)
            if v:
                p = resolve_path(v, root=root, data_dir=data_dir)
                self._optional[k] = p

        # hash inputs
        self._hash_inputs = [
            self._align_bam, self._bed_file, self._gtf_path,
            regions_tsv, align_pb.signature,
        ] + list(self._optional.values())
        return []

    def expected_outputs(self) -> Dict[str, Path]:
        return {
            "combined_bam":  Path("combined_region.bam"),
            "combined_bai":  Path("combined_region.bam.bai"),
            "combined_bed":  Path("combined_region.bed"),
            "combined_gtf":  Path("combined_region.gtf"),
            "slice_region_details": Path("slice_region_details.tsv"),
        }

    def _run(self):
        import time, platform
        from datetime import datetime, timezone

        sig = self.signature
        pb  = PathBuilder(self.work_dir, self.run_id, self.name, sig)
        d   = pb.stage_dir
        d.mkdir(parents=True, exist_ok=True)

        # BAM
        bam_out = d / "combined_region.bam"
        if Reinstate.decide(d, bam_out, needs_qc=False, stage_name=self.name) == "skip":
            return pb
        t0 = time.time()
        view = ["samtools","view","-b",str(self._align_bam),*self._intervals]
        sort = ["samtools","sort","-o",str(bam_out),"-T",str(d/"tmp_sort")]
        p_view = subprocess.Popen(view, stdout=subprocess.PIPE)
        try:
            subprocess.run(sort, stdin=p_view.stdout, check=True)
        finally:
            if p_view.stdout: p_view.stdout.close()
            p_view.wait()
        subprocess.run(["samtools","index",str(bam_out)], check=True)

        # BED
        bed_out = d / "combined_region.bed"
        if self._bed_file.exists():
            filter_file_by_regions(self._bed_file, bed_out, self._lookup, "bed")
        else:
            warnings.warn(f"BED missing: {self._bed_file}", UserWarning)

        # GTF
        gtf_out = d / "combined_region.gtf"
        if self._gtf_path.exists():
            filter_file_by_regions(self._gtf_path, gtf_out, self._lookup, "gtf")
        else:
            warnings.warn(f"GTF missing: {self._gtf_path}", UserWarning)

        # optional extras
        for k, pth in self._optional.items():
            out = d / f"combined_region.{k.replace('_','.')}"
            if pth.exists():
                # Determine filetype from suffix
                suf = pth.suffix.lower()
                if suf == ".bed":
                    filter_file_by_regions(pth, out, self._lookup, "bed")
                elif suf == ".tab":
                    filter_file_by_regions(pth, out, self._lookup, "tab")
                else:
                    warnings.warn(f"Skipping unsupported optional: {pth.name}", UserWarning)
            else:
                warnings.warn(f"Optional missing: {pth}", UserWarning)

        # slice_region_details.tsv
        det = d / "slice_region_details.tsv"
        with open(det, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["chrom","start","end","span_bp"])
            for c,s,e in self._regions:
                w.writerow([c,s,e, e-s+1])

        runtime = round(time.time()-t0,2)
        qc_metrics = {}
        if (fn := QC_REGISTRY.get(self.name)):
            qc_metrics = fn(det, d, runtime_sec=runtime)

        meta = {
            "stage": self.name,
            "signature": sig,
            "cmd": "(internal slice combined)",
            "exit_code": 0,
            "runtime_sec": runtime,
            "started": datetime.fromtimestamp(t0,timezone.utc).isoformat(),
            "ended": datetime.now(timezone.utc).isoformat(),
            "host": platform.node(),
            "qc": qc_metrics,
        }
        write_marker(pb, meta)
        return pb

    def run(self):
        self.build_cmd()
        self._input_hashes = hash_many(self._hash_inputs)
        self._flags_str = "combined-slice"
        return self._run()

