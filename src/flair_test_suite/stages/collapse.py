from __future__ import annotations

import logging
from pathlib import Path

from .base import StageBase
from .stage_utils import resolve_path, parse_cli_flags, get_stage_config
import logging
from ..lib.input_hash import hash_many


class CollapseStage(StageBase):
    """
    FLAIR collapse stage:
    - Consumes `slice` or `correct` stage output BED
    - Runs `flair collapse` with user-specified flags
    - Produces isoform BED/GTF (BED is primary output)
    """
    name = "collapse"
    requires: tuple[str, ...] = ()  # manual upstream handling
    primary_output_key = "isoforms_bed"

    def _locate_query_bed(self) -> tuple[Path, str, str]:
        # Prefer slice
        if "slice" in self.upstreams:
            slice_pb = self.upstreams["slice"]
            bed = slice_pb.stage_dir / "combined_region.bed"
            if bed.exists():
                return bed, slice_pb.signature, "slice"
        # Fallback to correct
        if "correct" in self.upstreams:
            corr_pb = self.upstreams["correct"]
            bed = corr_pb.stage_dir / f"{self.run_id}_all_corrected.bed"
            if bed.exists():
                return bed, corr_pb.signature, "correct"
        raise RuntimeError("collapse stage requires upstream `slice` or `correct` with an existing output BED.")

    def build_cmd(self) -> list[str]:
        cfg = self.cfg

        # Locate query BED, upstream signature, and branch
        query_bed, upstream_sig, branch = self._locate_query_bed()

        # --- resolve genome and reads ---
        if branch == "slice":
            prefix = query_bed.stem
            read_file = query_bed.parent / f"{prefix}.fa"
            if not read_file.exists():
                raise FileNotFoundError(f"Slice branch reads not found: {read_file}")
            raw_reads = read_file
        else:
            # correct branch: derive raw reads from config
            read_field = getattr(cfg.run, "reads_file", None) or getattr(cfg.run, "reads_files", None)
            if not read_field:
                raise RuntimeError("Config missing `reads_file(s)` required for collapse in correct branch.")
            raw_reads = read_field

        resolved = self.resolve_stage_inputs({
            "genome": cfg.run.genome_fa,
            "reads": raw_reads,
        })
        genome = resolved["genome"]
        reads = resolved["reads"]
        if isinstance(reads, (str, Path)):
            reads = [reads]
        self._genome_fa_abs = str(genome)  # For QC

        # --- parse flags and extra inputs ---
        flag_parts, extra_inputs = self.resolve_stage_flags()
        self._hash_inputs = [query_bed, genome, upstream_sig, *reads, *extra_inputs]
        self._flags_components = flag_parts

        if not flag_parts:
            logging.warning(
                "No extra flags configured for collapse; using defaults."
            )

        # --- build and return collapse command ---
        cmd = [
            "flair", "collapse",
            "-q", str(query_bed),
            "-g", str(genome),
            "-r", *map(str, reads),
            "-o", self.run_id,
            *flag_parts,
        ]
        logging.debug(f"[Collapse] Final command: {' '.join(map(str, cmd))}")
        return cmd

    def expected_outputs(self) -> dict[str, Path]:
        base = self.run_id
        return {
            "isoforms_bed": Path(f"{base}.isoforms.bed"),
            "isoforms_gtf": Path(f"{base}.isoforms.gtf"),
        }

    def collect_qc(self, pb):
        # Placeholder for future collapse QC metrics
        return {}
