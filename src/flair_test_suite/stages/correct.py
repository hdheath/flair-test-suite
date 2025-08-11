from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Tuple

import pandas as pd

from .base import StageBase
from .stage_utils import parse_cli_flags, get_stage_config
from ..lib.paths import PathBuilder
from ..lib.input_hash import hash_many
from ..qc.correct_qc import run_qc
from ..qc import write_metrics
from ..qc.qc_utils import count_lines


class CorrectStage(StageBase):
    """
    Runs `flair correct` on either:
      - the align BED (non-regionalized), or
      - each regionalized BED (regionalized),
    Skips empty inputs, and runs QC (regionalized = after-only motifs).
    """
    name = "correct"
    requires = ("align",)
    primary_output_key = "corrected"

    @property
    def tool_version(self) -> str:
        # Use configured run version for signature (or override to detect flair --version)
        return str(self.cfg.run.version)

    # -------- command builder (single source of truth) --------
    def build_cmds(self) -> List[List[str]]:
        cfg = self.cfg
        self._bed_files: List[Tuple[Path, str]] = []

        # Resolve genome path for later QC
        resolved = self.resolve_stage_inputs({
            "genome": cfg.run.genome_fa,
            "reads": getattr(cfg, "reads_file", None) or cfg.run.reads_file,
        })
        genome = resolved["genome"]
        self._genome_fa_abs = str(genome)

        # Discover upstreams
        upstream_sig = None
        if "regionalize" in self.upstreams:
            reg_pb = self.upstreams["regionalize"]
            upstream_sig = reg_pb.signature
            details_path = reg_pb.stage_dir / "region_details.tsv"
            if not details_path.exists():
                raise RuntimeError(f"Expected region_details.tsv not found: {details_path}")
            regions = pd.read_csv(details_path, sep="\t")
            for _, reg in regions.iterrows():
                chrom, start, end = reg["chrom"], reg["start"], reg["end"]
                bed_file = reg_pb.stage_dir / f"{chrom}_{start}_{end}.bed"
                if not bed_file.exists():
                    raise RuntimeError(f"Expected region output BED not found: {bed_file}")
                self._bed_files.append((bed_file, f"{chrom}_{start}_{end}"))
        elif "align" in self.upstreams:
            aln_pb = self.upstreams["align"]
            upstream_sig = aln_pb.signature
            bed_file = aln_pb.stage_dir / f"{self.run_id}_flair.bed"
            if not bed_file.exists():
                raise RuntimeError(f"Expected align output BED not found: {bed_file}")
            self._bed_files.append((bed_file, self.run_id))
        else:
            raise RuntimeError("No valid upstream found for correct stage.")

        # Parse flags, set signature inputs
        flag_parts, extra_inputs = self.resolve_stage_flags()
        self._flags_components = flag_parts
        self._hash_inputs = [genome] + [bf[0] for bf in self._bed_files] + extra_inputs
        if upstream_sig:
            self._hash_inputs.append(upstream_sig)

        # Capture align read count for metadata (used only in non-regionalized QC)
        align_meta = self.upstreams["align"].metadata
        self._n_input_reads = align_meta.get("n_input_reads", align_meta.get("n_total_reads"))

        # Compose commands, skipping empty inputs
        flair_version = str(cfg.run.version)
        major_version = int(flair_version.split(".")[0])

        cmds: List[List[str]] = []
        for bed_file, region_tag in self._bed_files:
            # hard skip empties (zero bytes OR zero records)
            if (not bed_file.exists()) or (bed_file.stat().st_size == 0) or (count_lines(bed_file) == 0):
                logging.warning(f"[correct] Skipping empty BED: {bed_file}")
                continue

            out_prefix = region_tag  # FLAIR writes {out_prefix}_all_corrected.bed
            cmd = ["flair", "correct", "-q", str(bed_file), "-o", out_prefix, *flag_parts]
            if major_version < 3:
                cmd.extend(["-g", str(genome)])
            cmds.append(cmd)

        return cmds

    # Optional legacy shim if something still calls build_cmd()
    def build_cmd(self) -> list[str]:
        cmds = self.build_cmds()
        return cmds[0] if cmds else []

    def expected_outputs(self) -> dict[str, Path]:
        # Pick a concrete primary: first region tag if present else run_id
        first_tag = (self._bed_files[0][1] if getattr(self, "_bed_files", []) else self.run_id)
        outputs = {
            "corrected":   Path(f"{first_tag}_all_corrected.bed"),
            "inconsistent": Path(f"{first_tag}_all_inconsistent.bed"),
            "qc_sidecar":  Path("correct_qc.tsv"),
        }
        # Optional per-region QC files (not used for reinstate, just helpful)
        if getattr(self, "_bed_files", []) and "regionalize" in self.upstreams:
            for _, tag in self._bed_files:
                outputs[f"{tag}_qc"] = Path(tag) / "correct_qc.tsv"
        return outputs

    # -------- QC orchestration (stage-specific) --------
    def _run_qc(self, stage_dir: Path, primary: Path, runtime: float | None) -> dict:
        """
        Run correct-stage QC:
          - regionalized: after-only motifs per region, no align comparison
          - non-regionalized: full before/after/diff
        Always writes an aggregate sidecar at stage root.
        """
        is_regionalized = "regionalize" in self.upstreams
        qc_metrics = run_qc(
            bed_files=self._bed_files,
            stage_dir=stage_dir,
            n_input_reads=(None if is_regionalized else self._n_input_reads),
            align_sig=(None if is_regionalized else self.upstreams["align"].signature),
            genome_fa=self._genome_fa_abs,
            runtime_sec=runtime,
            is_regionalized=is_regionalized,
        )

        # Aggregate sidecar (so Reinstate sees QC complete)
        total_qc_time = sum((m or {}).get("qc_runtime_sec", 0) for m in qc_metrics.values())
        write_metrics(stage_dir, "correct", {
            "regions": len(self._bed_files),
            "qc_runtime_sec": round(total_qc_time, 2),
        })
        return qc_metrics
