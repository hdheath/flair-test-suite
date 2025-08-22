from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from .base import StageBase
from .stage_utils import get_stage_config, resolve_path


@dataclass
class ManifestEntry:
    """Representation of one row in the FLAIR combine manifest."""

    sample: str
    isoform_type: str
    bed: Path
    fasta: Optional[Path] = None
    read_map: Optional[Path] = None


class CombineStage(StageBase):
    """FLAIR combine stage.

    Builds a multi-column manifest describing isoform files and runs
    ``flair combine`` to merge them. Collapse outputs are discovered
    automatically. Additional entries may be provided via the ``manifest``
    key under ``flags``. Each manifest row contains:

    ``sample``\t``type``\t``bed``\t``fasta``\t``read_map``

    where ``fasta`` and ``read_map`` are optional. Any extra keys in
    ``flags`` are treated as CLI options for ``flair combine``.
    """

    name = "combine"
    requires = ("collapse",)
    primary_output_key = "combined_bed"

    _manifest_rel = "manifest.tsv"
    _manifest_entries: List[ManifestEntry]

    def __init__(self, cfg, run_id: str, work_dir: Path, upstreams=None):
        super().__init__(cfg, run_id, work_dir, upstreams)
        self._manifest_entries = []

    @property
    def tool_version(self) -> str:  # pragma: no cover - simple accessor
        return str(self.cfg.run.version)

    def build_cmd(self) -> List[str]:
        cfg = self.cfg
        stage_cfg = get_stage_config(cfg, self.name)
        raw_flags = dict(getattr(stage_cfg, "flags", {}) or {})

        manifest_sources = raw_flags.pop("manifest", [])
        if isinstance(manifest_sources, (str, Path, dict)):
            manifest_sources = [manifest_sources]

        data_dir = Path(cfg.run.data_dir)

        entries: List[ManifestEntry] = []

        # Discover collapse outputs automatically
        collapse_pb = self.upstreams.get("collapse")
        if collapse_pb:
            collapse_dir = collapse_pb.stage_dir
            for bed in sorted(collapse_dir.glob("*.isoforms.bed")):
                sample = bed.name[:-len(".isoforms.bed")]
                fasta = bed.with_suffix(".fa")
                fasta = fasta if fasta.exists() else None
                read_map = bed.with_name(f"{sample}.isoform.read.map.txt")
                read_map = read_map if read_map.exists() else None
                entries.append(
                    ManifestEntry(sample=sample, isoform_type="isoform", bed=bed, fasta=fasta, read_map=read_map)
                )
            if not entries:
                logging.warning(f"[combine] collapse outputs missing in {collapse_dir}")

        # Append user-specified entries
        for item in manifest_sources:
            if isinstance(item, dict):
                bed = resolve_path(item.get("bed"), data_dir=data_dir)
                sample = item.get("sample")
                if not sample:
                    sample = bed.name
                    if sample.endswith(".isoforms.bed"):
                        sample = sample[:-len(".isoforms.bed")]
                    elif sample.endswith(".bed"):
                        sample = sample[:-len(".bed")]
                iso_type = item.get("type", "isoform")
                fasta = item.get("fasta")
                fasta = resolve_path(fasta, data_dir=data_dir) if fasta else None
                read_map = item.get("read_map")
                read_map = resolve_path(read_map, data_dir=data_dir) if read_map else None
                if not bed.exists():
                    logging.warning(f"[combine] manifest entry missing: {bed}")
                if fasta and not fasta.exists():
                    logging.warning(f"[combine] manifest fasta missing: {fasta}")
                if read_map and not read_map.exists():
                    logging.warning(f"[combine] manifest read_map missing: {read_map}")
                entries.append(ManifestEntry(sample, iso_type, bed, fasta, read_map))
            else:
                bed = resolve_path(item, data_dir=data_dir)
                sample = bed.name
                if sample.endswith(".isoforms.bed"):
                    sample = sample[:-len(".isoforms.bed")]
                elif sample.endswith(".bed"):
                    sample = sample[:-len(".bed")]
                iso_type = "isoform"
                fasta = bed.with_suffix(".fa")
                fasta = fasta if fasta.exists() else None
                read_map = bed.with_name(f"{sample}.isoform.read.map.txt")
                read_map = read_map if read_map.exists() else None
                if not bed.exists():
                    logging.warning(f"[combine] manifest entry missing: {bed}")
                entries.append(ManifestEntry(sample, iso_type, bed, fasta, read_map))

        # De-duplicate while preserving order by BED path
        seen: set[Path] = set()
        unique_entries: List[ManifestEntry] = []
        for e in entries:
            if e.bed not in seen:
                unique_entries.append(e)
                seen.add(e.bed)

        if not unique_entries:
            raise RuntimeError("[combine] no isoform files found for manifest")

        self._manifest_entries = unique_entries

        # parse remaining flags for CLI and signature
        flag_parts, extra_inputs = self.resolve_stage_flags(raw_flags)
        upstream_sigs = [collapse_pb.signature] if collapse_pb else []
        hash_inputs: List[Path] = []
        for e in unique_entries:
            hash_inputs.append(e.bed)
            if e.fasta:
                hash_inputs.append(e.fasta)
            if e.read_map:
                hash_inputs.append(e.read_map)
        self._hash_inputs = [*hash_inputs, *upstream_sigs, *extra_inputs]
        self._flags_components = flag_parts

        cmd = ["flair", "combine", "-i", self._manifest_rel, "-o", self.run_id]
        cmd.extend(flag_parts)
        return cmd

    def build_cmds(self) -> List[List[str]] | None:  # pragma: no cover - single command
        return None

    def expected_outputs(self) -> dict[str, Path]:
        base = f"{self.run_id}_combined"
        return {
            "combined_bed": Path(f"{base}.isoforms.bed"),
            "combined_gtf": Path(f"{base}.isoforms.gtf"),
        }

    def _prepare_stage_dir(self):  # pragma: no cover - exercised via tests
        pb, outputs, primary, needs_qc = super()._prepare_stage_dir()
        manifest_path = pb.stage_dir / self._manifest_rel
        with manifest_path.open("w") as fh:
            for e in self._manifest_entries:
                fasta = e.fasta or ""
                read_map = e.read_map or ""
                fh.write(f"{e.sample}\t{e.isoform_type}\t{e.bed}\t{fasta}\t{read_map}\n")
        return pb, outputs, primary, needs_qc
