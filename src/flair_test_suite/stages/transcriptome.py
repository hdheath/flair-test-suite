# src/flair_test_suite/stages/transcriptome.py
# --------------------------------------------
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Tuple

from .base import StageBase
from .stage_utils import read_region_details, make_flair_cmd

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
        """
        Returns: ([(bam_path, tag)], [upstream_signatures], mode)
        mode ∈ {"regionalized", "standard"}
        """
        pairs: List[Tuple[Path, str]] = []
        upstream_sigs: List[Path] = []

        # Prefer regionalize if present
        if "regionalize" in self.upstreams:
            mode = "regionalized"
            reg_pb = self.upstreams["regionalize"]
            upstream_sigs.append(reg_pb.signature)

            details = reg_pb.stage_dir / "region_details.tsv"
            if not details.exists():
                raise RuntimeError(
                    f"[transcriptome] region_details.tsv not found: {details}"
                )

            for chrom, start, end in read_region_details(details):
                tag = f"{chrom}_{start}_{end}"
                bam = reg_pb.stage_dir / f"{tag}.bam"
                if not bam.exists() or bam.stat().st_size == 0:
                    logging.warning(
                        f"[transcriptome] Missing/empty region BAM, skipping: {bam}"
                    )
                    continue
                pairs.append((bam, tag))

            if not pairs:
                raise RuntimeError("[transcriptome] No non-empty per-region BAMs discovered.")

            return pairs, upstream_sigs, mode

        # Fallback: align
        if "align" in self.upstreams:
            mode = "standard"
            aln_pb = self.upstreams["align"]
            upstream_sigs.append(aln_pb.signature)
            bam = aln_pb.stage_dir / f"{self.run_id}_flair.bam"
            if not bam.exists() or bam.stat().st_size == 0:
                raise RuntimeError(f"[transcriptome] Align BAM missing/empty: {bam}")
            pairs.append((bam, self.run_id))
            return pairs, upstream_sigs, mode

        raise RuntimeError("[transcriptome] requires upstream `regionalize` or `align`.")

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
        self._hash_inputs = [genome, *[bp[0] for bp in bam_pairs], *upstream_sigs, *extra_inputs]
        self._flags_components = flag_parts

        cmds: List[List[str]] = []

        if mode == "regionalized":
            # One transcriptome call per region
            for bam, tag in bam_pairs:
                cmds.append(
                    make_flair_cmd(
                        "transcriptome",
                        bam=bam,
                        genome=genome,
                        out=tag,
                        flags=flag_parts,
                    )
                )
                logging.info(f"[transcriptome] Scheduled per-region transcriptome: {tag}")
        else:
            # Single standard call
            bam, _ = bam_pairs[0]
            cmds.append(
                make_flair_cmd(
                    "transcriptome",
                    bam=bam,
                    genome=genome,
                    out=self.run_id,
                    flags=flag_parts,
                )
            )
            logging.info(f"[transcriptome] Scheduled standard transcriptome: {self.run_id}")

        logging.debug(f"[transcriptome] mode={mode} commands={len(cmds)}")
        return cmds

    # Legacy shim
    def build_cmd(self) -> list[str]:
        cmds = self.build_cmds()
        return cmds[-1] if cmds else []

    # ───────────────────── expected outputs ─────────────────────
    def expected_outputs(self) -> dict[str, Path]:
        """
        Provide a CONCRETE primary so Reinstate can test existence.
        - Regionalized: first region tag’s isoforms.bed
        - Standard: run_id.isoforms.bed
        """
        if "regionalize" in self.upstreams:
            first = self._first_tag or f"{self.run_id}"
            return {
                "isoforms_bed": Path(f"{first}.isoforms.bed"),
                "isoforms_gtf": Path(f"{first}.isoforms.gtf"),
                # (human hint only) pattern for the rest:
                "isoforms_bed_pattern": Path("{chrom}_{start}_{end}.isoforms.bed"),
            }
        else:
            base = self.run_id
            return {
                "isoforms_bed": Path(f"{base}.isoforms.bed"),
                "isoforms_gtf": Path(f"{base}.isoforms.gtf"),
            }

    def _run_qc(self, stage_dir: Path, primary: Path, runtime: float | None) -> dict:
        qc: dict = {}
        try:
            from ..qc.ted import collect as ted_collect
            ted_collect(stage_dir, self.cfg)
            qc["TED"] = {"tsv": str(stage_dir / "TED.tsv")}
        except Exception as e:  # pragma: no cover - logging only
            logging.warning(f"[transcriptome] TED QC failed: {e}")
        return qc


