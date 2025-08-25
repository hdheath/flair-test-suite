from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

from .base import StageBase
from .stage_utils import get_stage_config, resolve_path


class QuantifyStage(StageBase):
    """FLAIR quantify stage.

    Uses a user-inputted reads manifest to run ``flair quantify`` using isoform FASTA
    discovered from upstream ``combine``, ``transcriptome``, or ``collapse``
    outputs.
    """

    name = "quantify"
    requires: tuple[str, ...] = ()
    primary_output_key = "counts"

    _manifest_rel = "reads_manifest.tsv"

    _user_manifest: Optional[Path]
    _isoforms_fa: Optional[Path]

    def __init__(self, cfg, run_id: str, work_dir: Path, upstreams=None):
        super().__init__(cfg, run_id, work_dir, upstreams)
        self._user_manifest = None
        self._isoforms_fa = None

    @property
    def tool_version(self) -> str:  # pragma: no cover - simple accessor
        return str(self.cfg.run.version)

    # ------------------------------------------------------------------
    def _find_isoforms_fasta(self) -> Path:
        """Return isoform FASTA path from upstreams or raise if missing."""
        search = [
            ("combine", f"{self.run_id}.fa"),
            ("transcriptome", f"{self.run_id}.isoforms.fa"),
            ("collapse", f"{self.run_id}.isoforms.fa"),
        ]
        for stage, fname in search:
            pb = self.upstreams.get(stage)
            if not pb:
                continue
            cand = pb.stage_dir / fname
            if cand.exists():
                return cand
            # fallback: any .fa in upstream dir
            try:
                cand = next(pb.stage_dir.glob("*.fa"))
                if cand.exists():
                    return cand
            except StopIteration:
                pass
            logging.debug(f"[quantify] no isoforms FASTA found in {pb.stage_dir}")
        raise RuntimeError("[quantify] isoform FASTA not found in upstream outputs")

    # ------------------------------------------------------------------
    def build_cmd(self) -> List[str]:
        cfg = self.cfg
        stage_cfg = get_stage_config(cfg, self.name)
        raw_flags = dict(getattr(stage_cfg, "flags", {}) or {})

        manifest_src = raw_flags.pop("manifest", None)
        if manifest_src is None:
            raise RuntimeError("[quantify] manifest file is required; provide via flags.manifest")
        data_dir = Path(cfg.run.data_dir)

        isoforms_fa = self._find_isoforms_fasta()
        self._isoforms_fa = isoforms_fa

        seen: dict[str, tuple[str, str, Path]] = {}
        reads_paths: List[Path] = []

        self._user_manifest = resolve_path(manifest_src, data_dir=data_dir)
        if not self._user_manifest.exists():
            raise FileNotFoundError(f"[quantify] manifest file missing: {self._user_manifest}")
        for line in self._user_manifest.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                raise RuntimeError(f"[quantify] manifest line has fewer than 4 columns: {line}")
            sample, condition, batch, reads_str = parts[:4]
            reads_p = resolve_path(reads_str, data_dir=data_dir)
            if not reads_p.exists():
                raise FileNotFoundError(f"[quantify] reads file missing: {reads_p}")
            key = (condition, batch, reads_p)
            existing = seen.get(sample)
            if existing and existing != key:
                raise RuntimeError(f"[quantify] conflicting manifest entry for sample '{sample}'")
            seen[sample] = key
            reads_paths.append(reads_p)

        flag_parts, extra_inputs = self.resolve_stage_flags(raw_flags)

        upstream_sigs = []
        for name in ("combine", "transcriptome", "collapse"):
            pb = self.upstreams.get(name)
            if pb:
                upstream_sigs.append(pb.signature)

        hash_inputs: List[Path] = [isoforms_fa, *reads_paths, *extra_inputs, *upstream_sigs, self._user_manifest]
        self._hash_inputs = hash_inputs
        self._flags_components = flag_parts

        cmd = [
            "flair",
            "quantify",
            "-r",
            str(self._user_manifest),
            "-i",
            str(isoforms_fa),
            "-o",
            self.run_id,
        ]
        cmd.extend(flag_parts)
        return cmd

    def build_cmds(self) -> List[List[str]] | None:  # pragma: no cover - single command
        return None

    # ------------------------------------------------------------------
    def expected_outputs(self) -> dict[str, Path]:
        # Use the combined counts TSV as the primary artifact for quantify
        base = f"{self.run_id}"
        return {
            "counts": Path(f"{base}.counts.tsv"),
        }

    # ------------------------------------------------------------------
    def _prepare_stage_dir(self):  # pragma: no cover - exercised via tests
        pb, outputs, primary, needs_qc = super()._prepare_stage_dir()
        if not self._user_manifest:
            raise RuntimeError("[quantify] manifest file is required; provide via flags.manifest")
    # Use the user-supplied manifest in-place (do not copy into the stage dir).
    # The manifest must be accessible at the resolved path when the stage runs.
        return pb, outputs, primary, needs_qc
