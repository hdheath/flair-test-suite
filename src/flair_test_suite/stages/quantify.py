from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from .base import StageBase
from .stage_utils import get_stage_config, resolve_path


@dataclass
class ManifestEntry:
    """One row of the reads manifest."""

    sample: str
    condition: str
    batch: str
    reads: Path


class QuantifyStage(StageBase):
    """FLAIR quantify stage.

    Builds a reads manifest and runs ``flair quantify`` using isoform FASTA
    discovered from upstream ``combine``, ``transcriptome``, or ``collapse``
    outputs. If ``combine`` is upstream, a user-supplied manifest is required
    to describe reads for all samples used to generate the combined isoforms.
    """

    name = "quantify"
    requires: tuple[str, ...] = ()
    primary_output_key = "isoform_tpm"

    _manifest_rel = "reads_manifest.tsv"

    _append_entries: List[ManifestEntry]
    _user_manifest: Optional[Path]
    _isoforms_fa: Optional[Path]

    def __init__(self, cfg, run_id: str, work_dir: Path, upstreams=None):
        super().__init__(cfg, run_id, work_dir, upstreams)
        self._append_entries = []
        self._user_manifest = None
        self._isoforms_fa = None
        self._combine_upstream = "combine" in (upstreams or {})

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
        data_dir = Path(cfg.run.data_dir)

        isoforms_fa = self._find_isoforms_fasta()
        self._isoforms_fa = isoforms_fa

        # reset entries each build to avoid duplicates when build_cmd() is called multiple times
        self._append_entries = []
        seen: dict[str, tuple[str, str, Path]] = {}
        reads_paths: List[Path] = []

        # Parse user manifest if provided
        if manifest_src:
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
                if existing:
                    if existing != key:
                        raise RuntimeError(f"[quantify] conflicting manifest entry for sample '{sample}'")
                    continue  # identical duplicate
                seen[sample] = key
                reads_paths.append(reads_p)

        # Append run's reads
        resolved = self.resolve_stage_inputs({
            "reads": getattr(cfg.run, "reads_files", None) or getattr(cfg.run, "reads_file", None)
        })
        reads = resolved["reads"]
        if isinstance(reads, (str, Path)):
            reads = [reads]
        reads = [Path(r) for r in reads]
        for idx, r in enumerate(reads, start=1):
            sample = self.run_id if len(reads) == 1 else f"{self.run_id}_{idx}"
            condition = "condition1"
            batch = "batch1"
            key = (condition, batch, r)
            existing = seen.get(sample)
            if existing:
                if existing == key:
                    continue  # identical row
                raise RuntimeError(f"[quantify] sample '{sample}' already exists with different data")
            if not r.exists():
                raise FileNotFoundError(f"[quantify] reads file missing: {r}")
            seen[sample] = key
            reads_paths.append(r)
            self._append_entries.append(ManifestEntry(sample, condition, batch, r))

        flag_parts, extra_inputs = self.resolve_stage_flags(raw_flags)

        upstream_sigs = []
        for name in ("combine", "transcriptome", "collapse"):
            pb = self.upstreams.get(name)
            if pb:
                upstream_sigs.append(pb.signature)

        hash_inputs: List[Path] = [isoforms_fa, *reads_paths, *extra_inputs, *upstream_sigs]
        if self._user_manifest:
            hash_inputs.append(self._user_manifest)
        self._hash_inputs = hash_inputs
        self._flags_components = flag_parts

        cmd = [
            "flair",
            "quantify",
            "-r",
            self._manifest_rel,
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
        base = f"{self.run_id}"
        return {
            "isoform_tpm": Path(f"{base}.isoform.tpm.tsv"),
            "gene_counts": Path(f"{base}.gene.counts.tsv"),
            "isoform_counts": Path(f"{base}.isoform.counts.tsv"),
            "gene_tpm": Path(f"{base}.gene.tpm.tsv"),
        }

    # ------------------------------------------------------------------
    def _prepare_stage_dir(self):  # pragma: no cover - exercised via tests
        pb, outputs, primary, needs_qc = super()._prepare_stage_dir()
        if self._combine_upstream and not self._user_manifest:
            raise RuntimeError(
                "[quantify] combine upstream detected but no manifest provided; "
                "supply a manifest with reads for all samples"
            )
        manifest_path = pb.stage_dir / self._manifest_rel
        if self._user_manifest:
            content = self._user_manifest.read_text()
            manifest_path.write_text(content)
            if content and not content.endswith("\n"):
                with manifest_path.open("a") as fh:
                    fh.write("\n")
        with manifest_path.open("a" if self._user_manifest else "w") as fh:
            for e in self._append_entries:
                fh.write(f"{e.sample}\t{e.condition}\t{e.batch}\t{e.reads}\n")
        return pb, outputs, primary, needs_qc
