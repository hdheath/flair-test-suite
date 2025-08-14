from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Tuple

from .base import StageBase
from .stage_utils import read_region_details, make_flair_cmd
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

    def _resolve_bed_files(self) -> Tuple[List[Tuple[Path, str]], Path | None]:
        """Discover BED inputs from upstream stages."""

        bed_files: List[Tuple[Path, str]] = []
        upstream_sig: Path | None = None

        if "regionalize" in self.upstreams:
            reg_pb = self.upstreams["regionalize"]
            upstream_sig = reg_pb.signature
            details_path = reg_pb.stage_dir / "region_details.tsv"
            if not details_path.exists():
                raise RuntimeError(
                    f"Expected region_details.tsv not found: {details_path}"
                )
            for chrom, start, end in read_region_details(details_path):
                bed_file = reg_pb.stage_dir / f"{chrom}_{start}_{end}.bed"
                if not bed_file.exists():
                    raise RuntimeError(
                        f"Expected region output BED not found: {bed_file}"
                    )
                bed_files.append((bed_file, f"{chrom}_{start}_{end}"))
        elif "align" in self.upstreams:
            aln_pb = self.upstreams["align"]
            upstream_sig = aln_pb.signature
            bed_file = aln_pb.stage_dir / f"{self.run_id}_flair.bed"
            if not bed_file.exists():
                raise RuntimeError(f"Expected align output BED not found: {bed_file}")
            bed_files.append((bed_file, self.run_id))
        else:
            raise RuntimeError("No valid upstream found for correct stage.")

        return bed_files, upstream_sig

    # -------- command builder (single source of truth) --------
    def build_cmds(self) -> List[List[str]]:
        cfg = self.cfg
        self._bed_files, upstream_sig = self._resolve_bed_files()

        resolved = self.resolve_stage_inputs({
            "genome": cfg.run.genome_fa,
            "reads": getattr(cfg, "reads_file", None) or cfg.run.reads_file,
        })
        genome = resolved["genome"]
        self._genome_fa_abs = str(genome)

        flag_parts, extra_inputs = self.resolve_stage_flags()
        self._flags_components = flag_parts
        self._hash_inputs = [genome] + [bf[0] for bf in self._bed_files] + extra_inputs
        if upstream_sig:
            self._hash_inputs.append(upstream_sig)

        align_meta = self.upstreams["align"].metadata
        self._n_input_reads = align_meta.get("n_input_reads", align_meta.get("n_total_reads"))

        flair_version = str(cfg.run.version)
        major_version = int(flair_version.split(".")[0])

        cmds: List[List[str]] = []
        for bed_file, region_tag in self._bed_files:
            if (
                (not bed_file.exists())
                or (bed_file.stat().st_size == 0)
                or (count_lines(bed_file) == 0)
            ):
                logging.warning(f"[correct] Skipping empty BED: {bed_file}")
                continue

            out_prefix = region_tag
            cmd = make_flair_cmd(
                "correct",
                bed=bed_file,
                genome=(genome if major_version < 3 else None),
                out=out_prefix,
                flags=flag_parts,
            )
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
