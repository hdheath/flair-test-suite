# src/flair_test_suite/stages/collapse.py
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Tuple

from .base import StageBase
from ..qc.qc_utils import bed_is_empty
from .stage_utils import (
    collect_upstream_pairs,
    build_flair_cmds,
    isoform_expected_outputs,
    run_ted_qc,
    run_sqanti_qc,
)


logger = logging.getLogger(__name__)
# Force-load TED QC so Reinstate knows transcriptome has QC
try:
    from ..qc import ted as _force_import_ted  # noqa: F401
except Exception:
    pass


class CollapseStage(StageBase):
    """
    FLAIR collapse stage (per-region when regionalized, single otherwise).
    Produces per-region {chrom}_{start}_{end}.isoforms.{bed,gtf} when regionalized.
    """
    name = "collapse"
    requires: tuple[str, ...] = ()
    primary_output_key = "isoforms_bed"

    @property
    def tool_version(self) -> str:
        return str(self.cfg.run.version)

    # cache for expected_outputs / QC
    _mode: str | None = None
    _region_tags: List[str] = []
    _first_tag: str | None = None

    def _collect_corrected_inputs(self) -> Tuple[List[Tuple[Path, str]], List[Path], str]:
        return collect_upstream_pairs(
            "collapse",
            self.upstreams,
            self.run_id,
            "all_corrected.bed",
            "{chrom}_{start}_{end}_all_corrected.bed",
            file_check=lambda p: p.exists() and not bed_is_empty(p),
        )

    def build_cmds(self) -> List[List[str]]:
        cfg = self.cfg

        bed_pairs, upstream_sigs, mode = self._collect_corrected_inputs()
        self._mode = mode
        self._region_tags = [t for _, t in bed_pairs]
        self._first_tag = self._region_tags[0] if self._region_tags else None

        resolved = self.resolve_stage_inputs({
            "genome": cfg.run.genome_fa,
            "reads": getattr(cfg.run, "reads_files", None) or getattr(cfg.run, "reads_file", None),
        })
        genome = resolved["genome"]
        reads = resolved["reads"]
        if isinstance(reads, (str, Path)):
            reads = [reads]
        reads = [Path(r) for r in reads]

        flag_parts, extra_inputs = self.resolve_stage_flags()
        if "--generate_map" not in flag_parts:
            flag_parts.append("--generate_map")

        self._hash_inputs = [
            genome,
            *reads,
            *[bp[0] for bp in bed_pairs],
            *upstream_sigs,
            *extra_inputs,
        ]
        self._flags_components = flag_parts

        cmds = build_flair_cmds(
            "collapse",
            bed_pairs,
            genome,
            reads,
            self.run_id,
            flag_parts,
            regionalized=(mode == "regionalized"),
            use_bed=True,
        )

        logger.debug("mode=%s commands=%s", mode, len(cmds))
        return cmds

    def expected_outputs(self) -> dict[str, Path]:
        """
        Provide a CONCRETE primary so Reinstate can test existence.
        - Regionalized: first region tag’s isoforms.bed
        - Standard: run_id.isoforms.bed
        """
        return isoform_expected_outputs(
            self.run_id,
            self._first_tag,
            regionalized=(self._mode == "regionalized"),
        )

    def _run_qc(self, stage_dir: Path, primary: Path, runtime: float | None) -> dict:
        ted = run_ted_qc(self.name, stage_dir, self.cfg, self.upstreams)
        sqanti = run_sqanti_qc(self.name, stage_dir, self.cfg, self.upstreams)
        return {**ted, **sqanti}
