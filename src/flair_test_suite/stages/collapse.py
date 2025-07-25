from __future__ import annotations

import warnings
from pathlib import Path

from .base import StageBase
from .stage_utils import resolve_path, parse_cli_flags
import logging
from ..lib.input_hash import hash_many
from ..lib.signature import compute_signature


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

        # Resolve genome reference
        root = Path(cfg.run.input_root)
        data_dir = Path(cfg.run.data_dir)
        genome = resolve_path(cfg.run.genome_fa, root=root, data_dir=data_dir)
        if not genome.exists():
            warnings.warn(f"Genome FASTA not found: {genome}", UserWarning)

        # Determine reads based on branch
        if branch == "slice":
            # Expect a FASTA named <prefix>.fa next to combined_region.bed
            prefix = query_bed.stem
            read_file = query_bed.parent / f"{prefix}.fa"
            if not read_file.exists():
                raise FileNotFoundError(f"Slice branch reads not found: {read_file}")
            reads = [read_file]
        else:
            # correct branch: derive raw reads from config (same as align)
            read_field = getattr(cfg.run, "reads_file", None) or getattr(cfg.run, "reads_files", None)
            if not read_field:
                raise RuntimeError("Config missing `reads_file(s)` required for collapse in correct branch.")
            if isinstance(read_field, (str, Path)):
                raw_reads = [read_field]
            else:
                raw_reads = list(read_field)
            reads = [resolve_path(p, root=root, data_dir=data_dir) for p in raw_reads]
            for r in reads:
                if not r.exists():
                    warnings.warn(f"Read file not found: {r}", UserWarning)

        # Parse additional collapse flags
        raw_flags = next(st.flags for st in cfg.run.stages if st.name == "collapse")
        flag_parts, extra_inputs = parse_cli_flags(raw_flags, root=root, data_dir=data_dir)

        # Signature inputs: query_bed, genome, upstream_sig, reads, extra files
        self._hash_inputs = [query_bed, genome, upstream_sig, *reads, *extra_inputs]
        self._flags_components = flag_parts

        if not flag_parts:
            warnings.warn(
                "No extra flags configured for collapse; using defaults.",
                UserWarning
            )

        # Build and return collapse command
        env = cfg.run.conda_env
        cmd = [
            "conda", "run", "-n", env,
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
