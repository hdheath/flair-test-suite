from __future__ import annotations
from pathlib import Path
import csv, subprocess, logging
from collections import defaultdict
from typing import Dict, List, Tuple

from .base import StageBase
from ..lib import PathBuilder
from ..lib.input_hash import hash_many
from ..lib.signature import write_marker
from ..lib.reinstate import Reinstate
from .stage_utils import resolve_path, filter_file_by_regions
from ..qc import QC_REGISTRY

__all__ = ["SliceStage"]

Region = Tuple[str, int, int]


def _read_regions(tsv: Path) -> List[Region]:
    regions: List[Region] = []
    for ln, raw in enumerate(tsv.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            logging.warning(f"regions TSV line {ln}: fewer than 3 columns")
            continue
        try:
            chrom = parts[0]; s = int(parts[1]); e = int(parts[2])
        except ValueError:
            logging.warning(f"regions TSV line {ln}: non-integer coords")
            continue
        if s > e:
            s, e = e, s
        regions.append((chrom, s, e))
    if not regions:
        raise RuntimeError(f"No valid regions parsed from {tsv}")
    return regions


class IntervalLookup:
    def __init__(self, regions: List[Region]):
        idx: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
        for c, s, e in regions:
            idx[c].append((s, e))
        for c, lst in idx.items():
            lst.sort()
        self.idx = idx

    def contains(self, chrom: str, start: int, end: int) -> bool:
        if chrom not in self.idx:
            return False
        for s, e in self.idx[chrom]:
            if start >= s and end <= e:
                return True
            if start < s:
                break
        return False


class SliceStage(StageBase):
    """
    Slice: create region‑restricted BAM/BED/GTF/FASTA.

    Source BED priority:
      1. User override via [run.stages.flags].bed
      2. Corrected BED (if correct stage exists & file present)
      3. Align BED
    """
    name = "slice"
    requires = ("align", "correct")
    primary_output_key = "combined_bam"

    def build_cmd(self):
        cfg = self.cfg
        root = Path(cfg.run.input_root)
        data_dir = Path(cfg.run.data_dir)

        # Determine upstream BED/BAM
        align_pb = self.upstreams.get("align")
        correct_pb = self.upstreams.get("correct")

        # Set BAM and BED defaults
        self._align_bam = None
        self._align_bed = None

        if align_pb:
            self._align_bam = align_pb.stage_dir / f"{self.run_id}_flair.bam"
            self._align_bed = align_pb.stage_dir / f"{self.run_id}_flair.bed"
        if correct_pb:
            self._correct_bed = correct_pb.stage_dir / f"{self.run_id}_all_corrected.bed"

        flags = next(st.flags for st in cfg.run.stages if st.name == "slice")

        # GTF path (required)
        gtf = flags.get("gtf")
        if not gtf:
            raise RuntimeError("No GTF specified in slice stage flags.")
        self._gtf_path = resolve_path(gtf, root=root, data_dir=data_dir)

        # Optional user override
        override_bed = flags.get("bed")
        self._bed_source = "align"
        if override_bed:
            self._bed_file = resolve_path(override_bed, root=root, data_dir=data_dir)
            self._bed_source = "override"
        elif correct_pb and self._correct_bed.exists():
            self._bed_file = self._correct_bed
            self._bed_source = "corrected"
        elif self._align_bed and self._align_bed.exists():
            self._bed_file = self._align_bed
            self._bed_source = "align"
        else:
            raise RuntimeError("No valid BED file found for slice stage.")

        regions_tsv = flags.get("regions_tsv")
        if not regions_tsv:
            raise RuntimeError("No regions_tsv specified in slice stage flags.")
        regions_tsv_path = resolve_path(regions_tsv, root=root, data_dir=data_dir)
        self._regions = _read_regions(regions_tsv_path)
        self._intervals = [f"{c}:{s}-{e}" for c, s, e in self._regions]
        self._lookup = IntervalLookup(self._regions)

        # Optional extras
        opt_keys = [
            "junctions",
            "experiment_5_prime_regions_bed_file",
            "experiment_3_prime_regions_bed_file",
            "reference_5_prime_regions_bed_file",
            "reference_3_prime_regions_bed_file",
        ]
        self._optional: Dict[str, Path] = {}
        for k in opt_keys:
            v = flags.get(k)
            if v:
                self._optional[k] = resolve_path(v, root=root, data_dir=data_dir)

        # Hash inputs
        self._hash_inputs = [
            self._align_bam, self._bed_file, self._gtf_path,
            regions_tsv, align_pb.signature, *self._optional.values()
        ]
        return []

    def expected_outputs(self) -> dict[str, Path]:
        return {
            "combined_bam": Path("combined_region.bam"),
            "combined_bai": Path("combined_region.bam.bai"),
            "combined_bed": Path("combined_region.bed"),
            "combined_gtf": Path("combined_region.gtf"),
            "combined_fa": Path("combined_region.fa"),
            "slice_region_details": Path("slice_region_details.tsv"),
        }

    def _run(self):
        import time, platform
        from datetime import datetime, timezone

        sig = self.signature
        pb = PathBuilder(self.work_dir, self.run_id, self.name, sig)
        d = pb.stage_dir
        d.mkdir(parents=True, exist_ok=True)

        # Slice BAM
        bam_out = d / "combined_region.bam"
        if Reinstate.decide(d, bam_out, needs_qc=False, stage_name=self.name) == "skip":
            logging.info(f"[SKIP] {self.name} complete (sig={self.signature})")
            return pb
        t0 = time.time()
        view = ["samtools", "view", "-b", str(self._align_bam), *self._intervals]
        sort = ["samtools", "sort", "-o", str(bam_out), "-T", str(d / "tmp_sort")]
        p_view = subprocess.Popen(view, stdout=subprocess.PIPE)
        try:
            subprocess.run(sort, stdin=p_view.stdout, check=True)
        finally:
            if p_view.stdout: p_view.stdout.close()
            p_view.wait()
        subprocess.run(["samtools", "index", str(bam_out)], check=True)

        # FASTA
        with open(d / "combined_region.fa", "w") as fa:
            subprocess.run(["samtools", "fasta", str(bam_out)], stdout=fa, check=True)

        # BED slice
        bed_out = d / "combined_region.bed"
        if self._bed_file.exists():
            filter_file_by_regions(self._bed_file, bed_out, self._lookup, "bed")
        else:
            logging.warning(f"BED missing: {self._bed_file}")

        # GTF slice
        gtf_out = d / "combined_region.gtf"
        if self._gtf_path.exists():
            filter_file_by_regions(self._gtf_path, gtf_out, self._lookup, "gtf")
        else:
            logging.warning(f"GTF missing: {self._gtf_path}")

        # Extras
        for k, pth in self._optional.items():
            out = d / f"combined_region.{k.replace('_','.')}"
            if pth.exists():
                suf = pth.suffix.lower()
                if suf == ".bed":
                    filter_file_by_regions(pth, out, self._lookup, "bed")
                elif suf == ".tab":
                    filter_file_by_regions(pth, out, self._lookup, "tab")
                else:
                    logging.warning(f"Skipping unsupported optional: {pth.name}")
            else:
                logging.warning(f"Optional missing: {pth}")

        # Details
        det = d / "slice_region_details.tsv"
        with open(det, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["chrom", "start", "end", "span_bp"])
            for c, s, e in self._regions:
                w.writerow([c, s, e, e - s + 1])

        runtime = round(time.time() - t0, 2)
        qc_metrics = {}
        if (fn := QC_REGISTRY.get(self.name)):
            qc_metrics = fn(det, d, runtime_sec=runtime)

        meta = {
            "stage": self.name,
            "signature": sig,
            "bed_source": self._bed_source,
            "cmd": "(internal slice combined)",
            "exit_code": 0,
            "runtime_sec": runtime,
            "started": datetime.fromtimestamp(t0, timezone.utc).isoformat(),
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
