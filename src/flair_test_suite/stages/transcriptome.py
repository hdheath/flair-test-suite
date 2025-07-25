# src/flair_test_suite/stages/transcriptome.py
# --------------------------------------------
from __future__ import annotations

import warnings
from pathlib import Path

from .base import StageBase
from .stage_utils import resolve_path, parse_cli_flags


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

        # Genome
        root = Path(cfg.run.input_root)
        data_dir = Path(cfg.run.data_dir)
        genome = resolve_path(cfg.run.genome_fa, root=root, data_dir=data_dir)
        if not genome.exists():
            warnings.warn(f"Genome FASTA not found: {genome}", UserWarning)

        # Parse user flags (gtf, junction_tab, etc.)
        flags = next((st.flags for st in self.cfg.run.stages if st.name == self.name), None)
        raw_flags = flags or {}
        flag_parts, extra_inputs = parse_cli_flags(raw_flags, root=root, data_dir=data_dir)

        # Inject threads if missing
        threads = getattr(cfg.run, "threads", None)
        if threads and not any(f in ("-t", "--threads") for f in flag_parts):
            flag_parts.extend(["-t", str(threads)])

        # Signature inputs
        self._hash_inputs = [bam_path, genome, upstream_sig, *extra_inputs]
        self._flags_components = flag_parts

        if not flag_parts:
            warnings.warn(
                "No extra flags configured for transcriptome stage; using defaults.",
                UserWarning
            )

        env = cfg.run.conda_env
        cmd = [
            "conda", "run", "-n", env,
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

