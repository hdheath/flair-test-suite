from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Tuple

from .base import StageBase
from .stage_utils import collect_upstream_pairs, make_flair_cmd
from ..qc.correct_qc import run_qc
from ..qc import write_metrics
from ..qc.qc_utils import bed_is_empty


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

    def _resolve_bed_files(self) -> Tuple[List[Tuple[Path, str]], List[Path], str]:
        return collect_upstream_pairs(
            "correct",
            self.upstreams,
            self.run_id,
            "flair.bed",
            "{chrom}_{start}_{end}.bed",
        )

    # -------- command builder (single source of truth) --------
    def build_cmds(self) -> List[List[str]]:
        cfg = self.cfg
        self._bed_files, upstream_sigs, _ = self._resolve_bed_files()

        resolved = self.resolve_stage_inputs({
            "genome": cfg.run.genome_fa,
            "reads": getattr(cfg, "reads_file", None) or cfg.run.reads_file,
        })
        genome = resolved["genome"]
        self._genome_fa_abs = str(genome)

        flag_parts, extra_inputs = self.resolve_stage_flags()
        self._flags_components = flag_parts
        self._hash_inputs = [
            genome,
            *[bf[0] for bf in self._bed_files],
            *upstream_sigs,
            *extra_inputs,
        ]

        align_meta = self.upstreams["align"].metadata
        self._n_input_reads = align_meta.get("n_input_reads", align_meta.get("n_total_reads"))
        self._read_count_method = align_meta.get("read_count_method", "exact")

        flair_version = str(cfg.run.version)
        major_version = int(flair_version.split(".")[0])

        cmds: List[List[str]] = []
        for bed_file, region_tag in self._bed_files:
            if bed_is_empty(bed_file):
                logging.warning(f"[correct] Skipping missing/empty BED: {bed_file}")
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

    def expected_outputs(self) -> dict[str, Path]:
        # Pick a concrete primary: first region tag if present else run_id
        first_tag = (self._bed_files[0][1] if getattr(self, "_bed_files", []) else self.run_id)
        outputs = {
            "corrected":   Path(f"{first_tag}_all_corrected.bed"),
            "inconsistent": Path(f"{first_tag}_all_inconsistent.bed"),
            "qc_sidecar":  Path("qc/correct_qc.tsv"),
        }
        # Optional per-region QC files (not used for reinstate, just helpful)
        if getattr(self, "_bed_files", []) and "regionalize" in self.upstreams:
            for _, tag in self._bed_files:
                outputs[f"{tag}_qc"] = Path("qc") / tag / "correct_qc.tsv"
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
            read_count_method=(None if is_regionalized else self._read_count_method),
        )

        # Aggregate sidecar only for regionalized runs (per-region metrics already written)
        total_qc_time = sum((m or {}).get("qc_runtime_sec", 0) for m in qc_metrics.values())
        if is_regionalized:
            write_metrics(stage_dir, "correct", {
                "regions": len(self._bed_files),
                "qc_runtime_sec": round(total_qc_time, 2),
            })
        return qc_metrics
