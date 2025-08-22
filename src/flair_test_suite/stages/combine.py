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

    Builds a manifest describing isoform files and runs ``flair combine``
    to merge them. By default the collapse outputs from the current run are
    used. If a ``manifest`` flag is supplied it should point to an existing
    manifest TSV. That file is copied into the stage directory and the
    collapse outputs discovered in this run are appended to the end. Each
    manifest row contains ``sample``\t``type``\t``bed``\t``fasta``\t``read_map``
    where ``fasta`` and ``read_map`` are optional. Any extra keys in
    ``flags`` are treated as CLI options for ``flair combine``.
    """

    name = "combine"
    requires = ("collapse",)
    primary_output_key = "combined_bed"

    _manifest_rel = "manifest.tsv"
    _manifest_entries: List[ManifestEntry]
    _user_manifest: Optional[Path]

    def __init__(self, cfg, run_id: str, work_dir: Path, upstreams=None):
        super().__init__(cfg, run_id, work_dir, upstreams)
        self._manifest_entries = []
        self._user_manifest = None

    @property
    def tool_version(self) -> str:  # pragma: no cover - simple accessor
        return str(self.cfg.run.version)

    def build_cmd(self) -> List[str]:
        cfg = self.cfg
        stage_cfg = get_stage_config(cfg, self.name)
        raw_flags = dict(getattr(stage_cfg, "flags", {}) or {})

        manifest_src = raw_flags.pop("manifest", None)
        data_dir = Path(cfg.run.data_dir)

        entries: List[ManifestEntry] = []

        # Discover isoform outputs automatically from upstreams: prefer `collapse`,
        # but also include `transcriptome` outputs when present.
        collapse_pb = self.upstreams.get("collapse")
        transcriptome_pb = self.upstreams.get("transcriptome")

        def _collect_from(pb):
            collected = 0
            if not pb:
                return collected
            d = pb.stage_dir
            for bed in sorted(d.glob("*.isoforms.bed")):
                sample = bed.name[:-len(".isoforms.bed")]
                fasta = bed.with_suffix(".fa")
                fasta = fasta if fasta.exists() else None
                read_map = bed.with_name(f"{sample}.isoform.read.map.txt")
                read_map = read_map if read_map.exists() else None
                entries.append(
                    ManifestEntry(sample=sample, isoform_type="isoform", bed=bed, fasta=fasta, read_map=read_map)
                )
                collected += 1
            if collected == 0:
                logging.debug(f"[combine] no isoforms found in upstream dir {d}")
            return collected

        collected = 0
        collected += _collect_from(collapse_pb)
        collected += _collect_from(transcriptome_pb)
        if collected == 0:
            # only warn if neither upstream provided any entries; user manifest may still be used
            if collapse_pb or transcriptome_pb:
                logging.warning(f"[combine] no isoform outputs found in collapse/transcriptome upstreams")

        # Load user manifest file if provided
        if manifest_src:
            self._user_manifest = resolve_path(manifest_src, data_dir=data_dir)
            if not self._user_manifest.exists():
                logging.warning(f"[combine] manifest file missing: {self._user_manifest}")

        # Store collapse entries for later appending
        if not entries and not self._user_manifest:
            raise RuntimeError("[combine] no isoform files found for manifest")

        self._manifest_entries = entries

        # parse remaining flags for CLI and signature
        flag_parts, extra_inputs = self.resolve_stage_flags(raw_flags)
        upstream_sigs = []
        if collapse_pb:
            upstream_sigs.append(collapse_pb.signature)
        if transcriptome_pb:
            upstream_sigs.append(transcriptome_pb.signature)
        hash_inputs: List[Path] = []
        if self._user_manifest:
            hash_inputs.append(self._user_manifest)
        for e in entries:
            hash_inputs.append(e.bed)
            if e.fasta:
                hash_inputs.append(e.fasta)
            if e.read_map:
                hash_inputs.append(e.read_map)
        self._hash_inputs = [*hash_inputs, *upstream_sigs, *extra_inputs]
        self._flags_components = flag_parts

        cmd = ["flair", "combine", "--manifest", self._manifest_rel, "-o", self.run_id]
        cmd.extend(flag_parts)
        return cmd

    def build_cmds(self) -> List[List[str]] | None:  # pragma: no cover - single command
        return None

    def expected_outputs(self) -> dict[str, Path]:
        # FLAIR combine produces several artifacts: a merged BED, counts TSV,
        # a FASTA of merged isoforms, and an isoform map file. Use the
        # run_id as the base name to match the tool output.
        base = f"{self.run_id}"
        return {
            "combined_bed": Path(f"{base}.bed"),
            "combined_counts": Path(f"{base}.counts.tsv"),
            "combined_fasta": Path(f"{base}.fa"),
            "combined_map": Path(f"{base}.isoform.map.txt"),
        }

    def _prepare_stage_dir(self):  # pragma: no cover - exercised via tests
        pb, outputs, primary, needs_qc = super()._prepare_stage_dir()
        manifest_path = pb.stage_dir / self._manifest_rel
        seen: set[Path] = set()
        content = ""
        max_sample = 0
        if self._user_manifest:
            content = self._user_manifest.read_text()
            manifest_path.write_text(content)
            for line in content.splitlines():
                parts = line.split("\t")
                if len(parts) >= 3:
                    seen.add(Path(parts[2]))
                if parts and parts[0].isdigit():
                    max_sample = max(max_sample, int(parts[0]))
            mode = "a"
        else:
            mode = "w"

        with manifest_path.open(mode) as fh:
            if mode == "a" and content and not content.endswith("\n"):
                fh.write("\n")
            next_sample = max_sample + 1
            for e in self._manifest_entries:
                if e.bed in seen:
                    continue
                fasta = e.fasta or ""
                read_map = e.read_map or ""
                fh.write(
                    f"{next_sample}\t{e.isoform_type}\t{e.bed}\t{fasta}\t{read_map}\n"
                )
                next_sample += 1
        return pb, outputs, primary, needs_qc
