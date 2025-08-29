# src/flair_test_suite/stages/transcriptome.py
# --------------------------------------------
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Tuple

from .base import StageBase
from .stage_utils import (
    collect_upstream_pairs,
    build_flair_cmds,
    isoform_expected_outputs,
    run_ted_qc,
    run_sqanti_qc,
)


logger = logging.getLogger(__name__)

# Ensure TED QC is registered for this stage even if QC package isn't imported elsewhere
try:  # pragma: no cover - best effort, optional dependency
    from ..qc import ted as _force_import_ted  # noqa: F401
except Exception:  # pragma: no cover - missing heavy deps
    pass
class TranscriptomeStage(StageBase):
    """
    FLAIR transcriptome stage (v3+)

    Modes:
      • Regionalized path (regionalize → transcriptome):
          - Discover per-region BAMs: {chrom}_{start}_{end}.bam (from regionalize)
          - Run ONE `flair transcriptome` per region with -o {chrom}_{start}_{end}
          - Outputs: {chrom}_{start}_{end}.isoforms.{bed,gtf}

      • Standard path (align → transcriptome):
          - Use {run_id}_flair.bam from align
          - Single `flair transcriptome` with -o {run_id}
          - Outputs: {run_id}.isoforms.{bed,gtf}
    """
    name = "transcriptome"
    requires: tuple[str, ...] = ()  # dynamic upstreams
    primary_output_key = "isoforms_bed"

    # cache to provide a concrete primary and light metadata
    _mode: str | None = None
    _region_tags: List[str] = []
    _first_tag: str | None = None

    @property
    def tool_version(self) -> str:
        return str(self.cfg.run.version)

    # ───────────────────────── helpers ─────────────────────────
    def _collect_bams(self) -> Tuple[List[Tuple[Path, str]], List[Path], str]:
        return collect_upstream_pairs(
            "transcriptome",
            self.upstreams,
            self.run_id,
            "flair.bam",
            "{chrom}_{start}_{end}.bam",
        )

    # ─────────────────── command builder ───────────────────
    def build_cmds(self) -> List[List[str]]:
        cfg = self.cfg

        # Ensure FLAIR v3+
        flair_version = str(cfg.run.version)
        try:
            major = int(flair_version.split(".")[0])
        except Exception:
            major = 0
        if major < 3:
            raise RuntimeError(
                f"flair transcriptome requires FLAIR >= 3.0.0; configured '{flair_version}'."
            )

        # Discover inputs
        bam_pairs, upstream_sigs, mode = self._collect_bams()
        self._mode = mode
        self._region_tags = [t for _, t in bam_pairs]
        self._first_tag = self._region_tags[0] if self._region_tags else None

        # Resolve genome (reads not required for transcriptome)
        resolved = self.resolve_stage_inputs({
            "genome": cfg.run.genome_fa,
        })
        genome = resolved["genome"]
        self._genome_fa_abs = str(genome)

        # Parse flags
        flag_parts, extra_inputs = self.resolve_stage_flags()

        # Signature inputs: genome + all BAMs + upstream sigs + extra inputs
        self._hash_inputs = [
            genome,
            *[bp[0] for bp in bam_pairs],
            *upstream_sigs,
            *extra_inputs,
        ]
        self._flags_components = flag_parts

        cmds = build_flair_cmds(
            "transcriptome",
            bam_pairs,
            genome,
            None,
            self.run_id,
            flag_parts,
            regionalized=(mode == "regionalized"),
            use_bam=True,
        )

        logger.debug("mode=%s commands=%s", mode, len(cmds))
        return cmds

    # Legacy shim
    # ───────────────────── expected outputs ─────────────────────
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

    def expected_qc_files(self, stage_dir: Path) -> dict[str, Path] | None:
        # Expect TED TSV always; add SQANTI results when sqanti env configured
        paths: dict[str, Path] = {"TED": stage_dir / "qc" / "ted" / "TED.tsv"}
        # If regionalized, also require browser mapping and referenced files
        regional_beds = list(stage_dir.glob("*_*_*.isoforms.bed"))
        if regional_beds:
            region_map = stage_dir / "qc" / "ted" / "transcriptome_browser" / "region_map.json"
            paths["TED_browser_map"] = region_map
            try:
                import json
                if region_map.exists():
                    data = json.loads(region_map.read_text())
                    for tag, p in data.items():
                        paths[f"browser:{tag}"] = Path(p)
            except Exception:
                pass
        env = getattr(getattr(self.cfg, "run", object()), "sqanti_env", None)
        if env:
            paths["SQANTI_results"] = stage_dir / "qc" / "sqanti" / "sqanti_results.tsv"
        return paths

