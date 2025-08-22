from __future__ import annotations

import logging
from pathlib import Path
from typing import List

from .base import StageBase
from .stage_utils import get_stage_config, resolve_path


class CombineStage(StageBase):
    """FLAIR combine stage.

    Builds a manifest from isoform files and runs ``flair combine`` to merge
    them. Any paths listed in the ``manifest`` key of the stage's ``flags``
    block are added after the upstream collapse output, which is always
    included automatically when present. Additional keys in ``flags`` are
    treated as CLI options for ``flair combine``.
    """

    name = "combine"
    requires = ("collapse",)
    primary_output_key = "combined_bed"

    _manifest_rel = "manifest.tsv"
    _manifest_entries: List[Path]

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
        if isinstance(manifest_sources, (str, Path)):
            manifest_sources = [manifest_sources]

        data_dir = Path(cfg.run.data_dir)
        # Start with collapse output(s) if present
        manifest_paths: list[Path] = []
        collapse_pb = self.upstreams.get("collapse")
        if collapse_pb:
            collapse_dir = collapse_pb.stage_dir
            for p in sorted(collapse_dir.glob("*.isoforms.bed")):
                manifest_paths.append(p)
            if not manifest_paths:
                logging.warning(f"[combine] collapse outputs missing in {collapse_dir}")

        # Append user-specified entries
        for item in manifest_sources:
            p = resolve_path(item, data_dir=data_dir)
            manifest_paths.append(p)
            if not p.exists():
                logging.warning(f"[combine] manifest entry missing: {p}")

        # De-duplicate while preserving order
        seen: set[Path] = set()
        unique_manifest = []
        for p in manifest_paths:
            if p not in seen:
                unique_manifest.append(p)
                seen.add(p)

        if not unique_manifest:
            raise RuntimeError("[combine] no isoform files found for manifest")

        self._manifest_entries = unique_manifest

        # parse remaining flags for CLI and signature
        flag_parts, extra_inputs = self.resolve_stage_flags(raw_flags)
        upstream_sigs = [collapse_pb.signature] if collapse_pb else []
        self._hash_inputs = [*unique_manifest, *upstream_sigs, *extra_inputs]
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
            for p in self._manifest_entries:
                fh.write(f"{p}\n")
        return pb, outputs, primary, needs_qc
