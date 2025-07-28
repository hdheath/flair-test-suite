# src/flair_test_suite/stages/transcriptome.py
# --------------------------------------------
from __future__ import annotations

import logging
from pathlib import Path

from .base import StageBase
from .stage_utils import resolve_path, parse_cli_flags, get_stage_config


class TranscriptomeStage(StageBase):
    """
    FLAIR transcriptome stage (v3+):
    - Consumes the sliced BAM if slice stage present, else the align BAM
    - Runs `flair transcriptome` to detect isoforms
    - Produces isoform BED/GTF (BED is primary output)
    """
    name = "transcriptome"
    requires: tuple[str, ...] = ()  # manual upstream handling
    primary_output_key = "isoforms_bed"

    def _locate_bam(self) -> tuple[Path, str]:
        """
        Return (bam_path, upstream_signature).
        Priority: slice > align.
        """
        if "slice" in self.upstreams:
            slice_pb = self.upstreams["slice"]
            bam = slice_pb.stage_dir / "combined_region.bam"
            if bam.exists():
                return bam, slice_pb.signature

        if "align" in self.upstreams:
            align_pb = self.upstreams["align"]
            bam = align_pb.stage_dir / f"{self.run_id}_flair.bam"
            if bam.exists():
                return bam, align_pb.signature

        raise RuntimeError(
            "transcriptome stage requires upstream `slice` or `align` with an existing BAM."
        )

    def build_cmd(self) -> list[str]:
        cfg = self.cfg

        # Ensure version supports transcriptome
        flair_version = str(cfg.run.version)
        try:
            major_version = int(flair_version.split(".")[0])
        except Exception:
            major_version = 0
        if major_version < 3:
            raise RuntimeError(
                f"flair transcriptome requires FLAIR >= 3.0.0; version '{flair_version}' configured."
            )

        bam_path, upstream_sig = self._locate_bam()

        # --- resolve genome and reads ---
        raw_reads = getattr(cfg, "reads_file", None) or cfg.run.reads_file
        resolved = self.resolve_stage_inputs({
            "genome": cfg.run.genome_fa,
            "reads": raw_reads,
        })
        genome = resolved["genome"]
        reads = resolved["reads"]
        self._genome_fa_abs = str(genome)  # For QC

        # --- parse flags and extra inputs ---
        flag_parts, extra_inputs = self.resolve_stage_flags()
        self._hash_inputs = [bam_path, genome, upstream_sig, *extra_inputs]
        self._flags_components = flag_parts

        if not flag_parts:
            logging.warning(
                "No extra flags configured for transcriptome stage; using defaults."
            )

        cmd = [
            "flair", "transcriptome",
            "-b", str(bam_path),
            "-g", str(genome),
            "-o", self.run_id,
            *flag_parts,
        ]
        return cmd

    def expected_outputs(self) -> dict[str, Path]:
        base = self.run_id
        return {
            "isoforms_bed": Path(f"{base}.isoforms.bed"),
            "isoforms_gtf": Path(f"{base}.isoforms.gtf")
        }

    def collect_qc(self, pb):
        return {}

