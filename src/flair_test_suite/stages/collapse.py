# src/flair_test_suite/stages/collapse.py
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Tuple

from .base import StageBase
from ..qc.qc_utils import bed_is_empty
from .stage_utils import read_region_details, make_flair_cmd
...
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
        if "correct" not in self.upstreams:
            raise RuntimeError("collapse requires upstream `correct` stage.")

        corr_pb = self.upstreams["correct"]
        corr_dir = corr_pb.stage_dir
        upstream_sigs: List[Path] = [corr_pb.signature]

        pairs: List[Tuple[Path, str]] = []

        if "regionalize" in self.upstreams:
            mode = "regionalized"
            reg_pb = self.upstreams["regionalize"]
            upstream_sigs.append(reg_pb.signature)

            details = reg_pb.stage_dir / "region_details.tsv"
            if not details.exists():
                raise RuntimeError(f"Expected region_details.tsv not found: {details}")

            for chrom, start, end in read_region_details(details):
                tag = f"{chrom}_{start}_{end}"
                bed = corr_dir / f"{tag}_all_corrected.bed"
                if bed_is_empty(bed):
                    logging.warning(
                        f"[collapse] Missing or empty corrected BED, skipping: {bed}"
                    )
                    continue
                pairs.append((bed, tag))
        else:
            mode = "standard"
            bed = corr_dir / f"{self.run_id}_all_corrected.bed"
            if not bed.exists():
                raise RuntimeError(f"Expected corrected BED not found: {bed}")
            if bed_is_empty(bed):
                raise RuntimeError(f"Corrected BED is empty: {bed}")
            pairs.append((bed, self.run_id))

        if not pairs:
            raise RuntimeError("[collapse] No non-empty corrected BED inputs discovered.")

        return pairs, upstream_sigs, mode

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

        cmds: List[List[str]] = []

        if mode == "regionalized":
            for bed, tag in bed_pairs:
                cmds.append(
                    make_flair_cmd(
                        "collapse",
                        bed=bed,
                        genome=genome,
                        reads=reads,
                        out=tag,
                        flags=flag_parts,
                    )
                )
                logging.info(f"[collapse] Scheduled per-region collapse: {tag}")
        else:
            bed, _ = bed_pairs[0]
            cmds.append(
                make_flair_cmd(
                    "collapse",
                    bed=bed,
                    genome=genome,
                    reads=reads,
                    out=self.run_id,
                    flags=flag_parts,
                )
            )
            logging.info(f"[collapse] Scheduled standard collapse: {self.run_id}")

        logging.debug(f"[collapse] mode={mode} commands={len(cmds)}")
        return cmds

    def build_cmd(self) -> list[str]:
        cmds = self.build_cmds()
        return cmds[-1] if cmds else []

    def expected_outputs(self) -> dict[str, Path]:
        """
        Provide a CONCRETE primary so Reinstate can test existence.
        - Regionalized: first region tag’s isoforms.bed
        - Standard: run_id.isoforms.bed
        """
        if "regionalize" in self.upstreams:
            # If build_cmds hasn't run yet, fall back to a reasonable default
            first = self._first_tag or f"{self.run_id}"
            return {
                "isoforms_bed": Path(f"{first}.isoforms.bed"),
                "isoforms_gtf": Path(f"{first}.isoforms.gtf"),
                # (Optional) Pattern hint for humans; Reinstate shouldn't rely on this:
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

            ted_collect(stage_dir, self.cfg, upstreams=self.upstreams)
            qc["TED"] = {"tsv": str(stage_dir / "TED.tsv")}

        except Exception as e:
            logging.warning(f"[collapse] TED QC failed: {e}")
        return qc
