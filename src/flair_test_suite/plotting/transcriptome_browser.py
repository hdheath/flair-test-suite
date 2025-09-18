#!/usr/bin/env python3
"""
plot_genome_browser_v13.py (reads + isoform bars only)

- Keeps original row-packing logic:
  * per-isoform relative packing + global base scan with adjacency guards
  * reserved rows scaffold (label, collapsed_bar, reads_guard, tss/tts, spacer)
  * flips rows for display so row 0 is at the top
  * unassigned reads packed with 3-row guard and placed at the bottom
- Removes ALL GTF and ALL TSS/TTS code/plots from earlier versions
- Isoform bars: draw at reserved 'collapsed_bar' row (br)

JSON config keys (unchanged):
  bam, gtf, mapping, collapsed_isoforms, genome, outdir,
  gene_height, gene_row_height, read_row_height, fig_width, iso_kde
"""

import json
import argparse
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple
from collections import defaultdict
from itertools import chain

import pandas as pd
import pysam
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from intervaltree import IntervalTree
from matplotlib.collections import PatchCollection, LineCollection


@dataclass(frozen=True)
class Config:
    gtf: Path
    bam: Optional[Path] = None
    genome: str = ""
    outdir: Path = Path(".")
    mapping: Optional[Path] = None
    collapsed_isoforms: Optional[Path] = None
    gene_height: float = 1.5
    gene_row_height: float = 3.0      # kept for backward compat; not used now
    read_row_height: float = 0.02
    fig_width: float = 12.0
    iso_kde: bool = True              # kept for compat; not used now
    max_depth_capacity: int = 1000     # cap figure height scaling to this many rows

    @staticmethod
    def from_json(p: Path) -> "Config":
        cfg = json.loads(p.read_text())
        return Config(
            bam=Path(cfg["bam"]) if cfg.get("bam") else None,
            gtf=Path(cfg["gtf"]),
            genome=cfg.get("genome", ""),
            outdir=Path(cfg.get("outdir", ".")),
            mapping=Path(cfg["mapping"]) if cfg.get("mapping") else None,
            collapsed_isoforms=Path(cfg["collapsed_isoforms"]) if cfg.get("collapsed_isoforms") else None,
            gene_height=cfg.get("gene_height", 1.5),
            gene_row_height=cfg.get("gene_row_height", 3.0),
            read_row_height=cfg.get("read_row_height", 0.02),
            fig_width=cfg.get("fig_width", 12.0),
            iso_kde=cfg.get("iso_kde", True),
            max_depth_capacity=int(cfg.get("max_depth_capacity", 1000)),
        )


def parse_args():
    p = argparse.ArgumentParser(
        description="Reads + collapsed isoform bars (no genes, no TSS/TTS).",
    )
    p.add_argument(
        "-c", "--config", type=Path, required=True,
        help=("JSON config keys: (bam), gtf, mapping, collapsed_isoforms, genome, outdir, "
              "gene_height, gene_row_height, read_row_height, fig_width, iso_kde"),
    )
    p.add_argument("--region", help="Limit plotting region (e.g. chr1:100-200)")
    return p.parse_args()


def pleasant_colors(seq):
    def ok(rgb):
        mx, mn = max(rgb), min(rgb)
        return (mx < 0.9) and not (mx - mn < 0.2 and mx > 0.7)
    return [c for c in seq if ok(c)]


def read_strand(rd: pysam.AlignedSegment) -> str:
    xs = rd.get_tag("XS") if rd.has_tag("XS") else None
    if xs in ("+", "-"):
        return xs
    return "-" if rd.is_reverse else "+"


def overlaps(tree: IntervalTree, blks) -> bool:
    return any(tree.overlap(s, e) for s, e in blks)


def add_blocks(tree: IntervalTree, blks) -> None:
    for s, e in blks:
        tree.addi(s, e, True)


def ensure_rows(occupancy, upto):
    missing = upto - (len(occupancy) - 1)
    for _ in range(max(0, missing)):
        occupancy.append(IntervalTree())


def reserve_iso_rows(reads_base: int, max_rel: int) -> dict[str, int]:
    """Reserve rows for one isoform stack.

    The ``reads_base`` is the first row index where reads are placed.
    We place the label and collapsed bar ABOVE the reads (smaller indices),
    and TSS/TTS/spacer BELOW the reads (larger indices).
    Visual order (top→bottom, after flip):
      label → collapsed_bar → reads → tss_panel → tts_panel → spacer
    """
    return {
        # Above the read stack
        "label":         reads_base - 4,
        "collapsed_bar": reads_base - 1,
        # Guard row below the read stack
        "reads":         reads_base + max_rel + 1,
        # Panels and spacer further below
        "tss_panel":     reads_base + max_rel + 6,
        "tts_panel":     reads_base + max_rel + 9,
        "spacer":        reads_base + max_rel + 12,
    }



def load_mapping(map_path: Optional[Path]):
    m: dict[str, str] = {}
    if not map_path or not map_path.exists():
        logging.warning(f"Mapping not found: {map_path}")
        return m
    with map_path.open() as f:
        for line in f:
            iso, reads = line.strip().split('\t')
            for r in reads.split(','):
                r = r.strip()
                if r:
                    m[r] = iso
    if not m:
        logging.warning("Mapping is empty.")
    return m


def assign_read_rows(blocks_list):
    """
    Row assignment using interval trees for fast overlap checks.
    Returns a list of row indices and total number of rows.
    """
    rows = []  # list of IntervalTrees
    assign = []
    for blks in blocks_list:
        placed = False
        for ridx, tree in enumerate(rows):
            # ensure none of these blocks overlap existing intervals
            if all(len(tree.overlap(b[0], b[1])) == 0 for b in blks):
                # insert blocks into this row's tree
                for b in blks:
                    tree[b[0]:b[1]] = True
                assign.append(ridx)
                placed = True
                break
        if not placed:
            tree = IntervalTree()
            for b in blks:
                tree[b[0]:b[1]] = True
            rows.append(tree)
            assign.append(len(rows) - 1)
    return assign, len(rows)


def _parse_region(region: Optional[str]) -> Tuple[Optional[Tuple[str, int, int]], bool]:
    """
    Parse a region string like ``"chr1:100-200"`` into a tuple and flag
    whether the span is too large for plotting.

    Returns a tuple ``(region_tuple, skip_plot)`` where ``region_tuple`` is
    ``(chrom, start, end)`` if parsing succeeds or ``None`` otherwise, and
    ``skip_plot`` is ``True`` when the span exceeds the 20kb plotting limit.
    """
    if not region:
        return None, False
    m = re.match(r"^([^:]+):(\d+)-(\d+)$", region.replace(",", ""))
    if not m:
        logging.warning(f"Could not parse region string: {region}")
        return None, False
    chrom, start_s, end_s = m.groups()
    start_i, end_i = int(start_s), int(end_s)
    # Do not decide to skip plotting here; leave span-based gating to the
    # caller (e.g., TED) so that policy is centralized. Always return the
    # parsed region tuple and let the caller decide whether to invoke the
    # plotting routine for very large regions.
    return (chrom, start_i, end_i), False


def generate(cfg: Config, region: Optional[str] = None) -> Optional[Path]:
    """Generate plot from configuration (reads + isoform bars only).

    Returns the Path to the saved PNG on success, or ``None`` if the plot was
    skipped (e.g., due to missing inputs).
    """
    # ---- Parse and validate region ----
    region_tuple, skip_plot = _parse_region(region)
    if region and region_tuple is None:
        return
    chrom: Optional[str] = None
    r0 = r1 = None
    if region_tuple:
        chrom, r0, r1 = region_tuple

    # ---- Resolve inputs (BAM, FASTA, mapping, isoform BED) ----
    def _resolve_bam_and_fa(gtf: Path, bam_cfg: Optional[Path], region_tuple: Optional[Tuple[str, int, int]]):
        bam: Optional[Path] = None
        reg_fa: Optional[Path] = None
        attempts: list[Path] = []
        if region_tuple:
            c, s, e = region_tuple
            tag = f"{c}_{s}_{e}"
            run_root = gtf.parent.parent.parent
            reg_root = run_root / "region_test"
            if reg_root.exists():
                for d in reg_root.iterdir():
                    cand = d / f"{tag}.bam"
                    if cand.exists():
                        bam = cand
                        fa_cand = d / f"{tag}.fa"
                        attempts.append(fa_cand)
                        if fa_cand.exists():
                            reg_fa = fa_cand
                        break
        if bam is None:
            bam = bam_cfg
        if region_tuple and (reg_fa is None) and bam is not None:
            try:
                c, s, e = region_tuple
                cand_fa = bam.parent / f"{c}_{s}_{e}.fa"
                attempts.append(cand_fa)
                if cand_fa.exists():
                    reg_fa = cand_fa
            except Exception:
                pass
        return bam, reg_fa, attempts

    gtf = cfg.gtf
    bam, reg_fa, _fa_attempts = _resolve_bam_and_fa(gtf, cfg.bam, region_tuple)

    if not bam or not bam.exists():
        logging.warning(f"BAM not found: {bam}")
        return
    if not gtf.exists():
        logging.warning(f"GTF not found (only used for regionalized BAM lookup): {gtf}")

    mapping_file = cfg.mapping if cfg.mapping and cfg.mapping.exists() else None
    if cfg.mapping and not mapping_file:
        logging.warning(f"Mapping not found: {cfg.mapping}")
    cis_bed = cfg.collapsed_isoforms if cfg.collapsed_isoforms and cfg.collapsed_isoforms.exists() else None
    if cfg.collapsed_isoforms and not cis_bed:
        logging.warning(f"Collapsed isoforms BED not found: {cfg.collapsed_isoforms}")

    if region_tuple:
        if reg_fa and reg_fa.exists():
            logging.info(f"[browser] Using regional FASTA: {reg_fa}")
        else:
            tried = ", ".join(str(p) for p in _fa_attempts) if _fa_attempts else "<none>"
            logging.warning(f"[browser] FASTA sequence not found; searched: {tried}")

    genome, outdir = cfg.genome, cfg.outdir
    gene_h, read_row_h, fig_w = cfg.gene_height, cfg.read_row_height, cfg.fig_width
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Load reads ----
    def _load_reads(bam: Path, region_tuple: Optional[Tuple[str, int, int]]):
        reads_blocks, reads_introns, read_strands, read_names = [], [], [], []
        contig_for_title = None
        r0_local = r0; r1_local = r1
        with pysam.AlignmentFile(bam, "rb") as bf:
            iterator = bf.fetch(*region_tuple) if region_tuple else bf.fetch()
            for rd in iterator:
                if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                    continue
                if contig_for_title is None:
                    contig_for_title = bf.get_reference_name(rd.reference_id)
                blks = rd.get_blocks()
                if not blks:
                    continue
                if r0_local is not None:
                    clipped = []
                    for s, e in blks:
                        s2, e2 = max(s, r0_local), min(e, r1_local)
                        if e2 > s2:
                            clipped.append((s2, e2))
                    blks = clipped
                    if not blks:
                        continue
                reads_blocks.append(blks)
                reads_introns.append([(blks[i][1], blks[i + 1][0]) for i in range(len(blks) - 1)])
                read_strands.append(read_strand(rd))
                read_names.append(rd.query_name)
        return reads_blocks, reads_introns, read_strands, read_names, contig_for_title

    reads_blocks, reads_introns, read_strands, read_names, contig_for_title = _load_reads(bam, region_tuple)
    if not reads_blocks:
        logging.warning("No reads loaded.")
        return
    all_starts = [b[0] for blks in reads_blocks for b in blks]
    all_ends = [b[1] for blks in reads_blocks for b in blks]
    read_min, read_max = min(all_starts), max(all_ends)

    # ---- Mapping and collapsed isoforms ----
    mapping = load_mapping(mapping_file)

    def _load_collapsed(cis_bed: Optional[Path]):
        cis_df = pd.DataFrame(columns=[
            "chr","iso_id","start","end","strand","blockCount","blockSizes","blockStarts"
        ])
        iso_blocks: dict[str, list[tuple[int, int]]] = {}
        if cis_bed and cis_bed.exists():
            tmp = pd.read_csv(cis_bed, sep="\t", header=None, comment="#")
            if tmp.shape[1] < 12:
                logging.warning("Collapsed isoforms BED missing BED12 fields.")
                tmp = tmp.iloc[:, :min(tmp.shape[1], 6)]
                cols = ["chr","start","end","iso_id","score","strand"][:tmp.shape[1]]
                tmp.columns = cols
                tmp["blockCount"], tmp["blockSizes"], tmp["blockStarts"] = np.nan, "", ""
            else:
                tmp = tmp.iloc[:, :12]
                tmp.columns = [
                    "chr","start","end","iso_id","score","strand",
                    "thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"
                ]
            cis_df = tmp[["chr","iso_id","start","end","strand","blockCount","blockSizes","blockStarts"]]
            for _, r in cis_df.iterrows():
                sizes = [int(x) for x in str(r.blockSizes).rstrip(",").split(",") if x]
                starts = [int(x) for x in str(r.blockStarts).rstrip(",").split(",") if x]
                if sizes and starts:
                    blks = [(int(r.start) + st, int(r.start) + st + sz) for st, sz in zip(starts, sizes)]
                else:
                    blks = [(int(r.start), int(r.end))]
                iso_blocks[r.iso_id] = blks
        return cis_df, iso_blocks

    cis_df, iso_blocks = _load_collapsed(cis_bed)

    # Group reads by isoform id (from mapping)
    iso_to_idxs, unassigned = defaultdict(list), []
    for i, name in enumerate(read_names):
        iso = mapping.get(name)
        (iso_to_idxs[iso].append(i) if iso else unassigned.append(i))

    base_colors = list(chain(plt.cm.tab20.colors, plt.cm.Set3.colors, plt.cm.Dark2.colors))
    pal = pleasant_colors(base_colors)
    isos = sorted(iso_to_idxs, key=lambda i: len(iso_to_idxs[i]), reverse=True)
    iso_colors = {iso: pal[i % len(pal)] for i, iso in enumerate(isos)}
    ua_color = 'lightgrey'

    # ---- Layout rows (packing) ----
    def _layout(reads_blocks, iso_to_idxs, cis_df):
        occupancy: list[IntervalTree] = []
        row_assign: list[Optional[int]] = [None] * len(reads_blocks)
        collapsed: dict[str, tuple[int, int, int, int, int, int]] = {}

        ih_local = gene_h * 0.6  # keep computation here to mirror original
        _ = ih_local  # silence linter; returned layout does not need this variable

        READS_OFFSET = 3  # leave room above reads for label + bar

        for iso in isos:
            idxs = iso_to_idxs[iso]
            sorted_idxs = sorted(idxs, key=lambda i: reads_blocks[i][0][0])
            blks_list, _ = zip(*[(reads_blocks[i], None) for i in sorted_idxs])
            rel, _ = assign_read_rows(blks_list)

            st = en = None
            if not cis_df.empty and iso in set(cis_df["iso_id"]):
                row = cis_df.loc[cis_df["iso_id"] == iso].iloc[0]
                st, en = int(row.start), int(row.end)
            if st is None:
                logging.warning(f"No collapsed bar for {iso}")
                unassigned.extend(idxs)
                continue

            left = min(reads_blocks[j][0][0] for j in idxs)
            right = max(reads_blocks[j][-1][1] for j in idxs)
            max_rel = max(rel)

            for base in range(len(occupancy) + 1):
                ok = True
                reads_base = base + READS_OFFSET
                for ri, blks in zip(rel, blks_list):
                    for off in (0, 1):
                        tgt = reads_base + ri + off
                        if tgt < len(occupancy) and overlaps(occupancy[tgt], blks):
                            ok = False
                            break
                    if not ok:
                        break
                if not ok:
                    continue
                spec_rows = reserve_iso_rows(reads_base, max_rel)
                for key, pr in spec_rows.items():
                    # Reserve full isoform span for all reserved rows to avoid
                    # downstream visual overlaps (labels/KDE/panels). This is a
                    # simplification over the prior thin-slice approach.
                    if key in ("label", "tss_panel", "tts_panel", "spacer", "collapsed_bar"):
                        seg = (st, en)
                    else:  # reads guard reserves full read span
                        seg = (left, right)
                    # Check this row and immediate neighbors to maintain a
                    # vertical guard band so thick rectangles do not collide
                    # with adjacent features on neighboring rows.
                    rows_to_check = [pr]
                    if key in ("label", "tss_panel", "tts_panel", "spacer", "collapsed_bar"):
                        rows_to_check += [pr - 1, pr + 1]
                    for rr in rows_to_check:
                        if rr < 0:
                            continue
                        if rr < len(occupancy) and overlaps(occupancy[rr], [seg]):
                            ok = False
                            break
                    if not ok:
                        break
                if ok:
                    break
            else:
                base = len(occupancy)
                reads_base = base + READS_OFFSET
                spec_rows = reserve_iso_rows(reads_base, max_rel)

            # Ensure room for all reserved rows plus a one-row guard below
            ensure_rows(occupancy, spec_rows["spacer"] + 1)
            for ri, idx in zip(rel, sorted_idxs):
                r = reads_base + ri
                row_assign[idx] = r
                add_blocks(occupancy[r], reads_blocks[idx])
            for key in ("label", "collapsed_bar", "reads", "tss_panel", "tts_panel", "spacer"):
                r = spec_rows[key]
                # Reserve full span for label/panels/bar to prevent overlap.
                if key in ("label", "tss_panel", "tts_panel", "spacer", "collapsed_bar"):
                    seg = (st, en)
                else:  # reads guard row blocks full read span
                    seg = (left, right)
                add_blocks(occupancy[r], [seg])
                # Add a guard band above/below for thick features
                if key in ("label", "tss_panel", "tts_panel", "spacer", "collapsed_bar"):
                    if r - 1 >= 0:
                        add_blocks(occupancy[r - 1], [seg])
                    ensure_rows(occupancy, r + 1)
                    add_blocks(occupancy[r + 1], [seg])

            collapsed[iso] = (
                spec_rows["label"],
                spec_rows["collapsed_bar"],
                spec_rows["reads"],
                spec_rows["tss_panel"],
                spec_rows["tts_panel"],
                spec_rows["spacer"],
            )

        # Unassigned reads (3-row guard), always placed in a dedicated
        # bottom block appended AFTER all isoform rows to avoid scattering.
        ua_start = len(occupancy)
        for idx in unassigned:
            blks = reads_blocks[idx]
            r = ua_start
            while True:
                # Ensure candidate row and a one-row neighbor below exist
                ensure_rows(occupancy, r + 1)
                # 3-row guard: check r-1, r, r+1
                bad = False
                # Above neighbor (only if exists)
                if r - 1 >= 0 and overlaps(occupancy[r - 1], blks):
                    bad = True
                # Current row
                if not bad and overlaps(occupancy[r], blks):
                    bad = True
                # Below neighbor
                if not bad and overlaps(occupancy[r + 1], blks):
                    bad = True
                if bad:
                    r += 1
                    continue
                # Place here
                row_assign[idx] = r
                add_blocks(occupancy[r], blks)
                break

        tot = len(occupancy)
        return row_assign, collapsed, tot

    row_assign, collapsed, tot = _layout(reads_blocks, iso_to_idxs, cis_df)

    # BED spans per iso id (for drawing model when blocks missing)
    ci_map = {r.iso_id: (int(r.start), int(r.end)) for _, r in cis_df.iterrows()} if not cis_df.empty else {}

    # ---- FASTA sequence (optional; for internal priming highlighting) ----
    def _load_fasta_seq(fa_path: Optional[Path]) -> str:
        if not fa_path or not fa_path.exists():
            return ""
        parts = []
        with fa_path.open() as fh:
            for ln in fh:
                if ln.startswith('>'):
                    continue
                parts.append(ln.strip())
        return ''.join(parts).upper()

    seq = _load_fasta_seq(reg_fa)

    # ---- Plot (reads + isoform bars) ----
    # Hard-code figure size to golden ratio with fixed height = 8 inches
    PHI = (1 + 5 ** 0.5) / 2
    fig_h = 8.0
    fig_w = fig_h * PHI
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=1200)
    ax = fig.add_subplot(1, 1, 1)
    for side in ['left', 'right', 'bottom', 'top']:
        ax.spines[side].set_visible(False)

    # Hard-coded compact heights for consistent spacing
    LABEL_H  = 0.8   # label rectangle total height
    BAR_H    = 0.5   # isoform bar height
    READ_H   = 0.5   # read block height (match reference style)
    INTRON_H = 0.05  # thin baseline/intron connector height
    KDE_H    = 1.0   # per-iso KDE panel height

    ih, rh = BAR_H, READ_H

    # Helper to flip raw row indices for display so row 0 is at the top.
    def disp_row(r: int) -> int:
        return (tot - 1 - r)

    # 1) Isoform labels + bars
    for iso in isos:
        if iso not in collapsed:
            continue
        label_r, br, reads_guard, tss_r, tts_r, spacer_r = collapsed[iso]
        if iso in ci_map:
            cs, ce = ci_map[iso]
        else:
            idxs = iso_to_idxs[iso]
            cs = min(reads_blocks[j][0][0] for j in idxs)
            ce = max(reads_blocks[j][-1][1] for j in idxs)
        if r0 is not None:
            cs = max(cs, r0); ce = min(ce, r1)
            if ce <= cs:
                continue
        blks_iso = iso_blocks.get(iso, [(cs, ce)])
        if r0 is not None:
            blks_iso = [(max(bs, r0), min(be, r1)) for (bs, be) in blks_iso if be > r0 and bs < r1]
            if not blks_iso:
                continue
        color = iso_colors.get(iso, 'black')
        bar_y, label_y = disp_row(br), disp_row(label_r)
        span_left = min(bs for bs, _ in blks_iso); span_right = max(be for _, be in blks_iso)
        for bs, be in blks_iso:
            ax.add_patch(Rectangle((bs, bar_y - ih/2), be - bs, ih,
                                   facecolor=color, edgecolor='black', linewidth=0.1, zorder=3))
        intr = [(blks_iso[i][1], blks_iso[i+1][0]) for i in range(len(blks_iso) - 1)]
        if intr:
            ax.add_collection(LineCollection([((a, bar_y), (b, bar_y)) for a, b in intr],
                                             colors='black', linewidths=0.5, zorder=3))
        iso_label = iso
        ax.add_patch(Rectangle((span_left, label_y - LABEL_H/2), span_right - span_left, LABEL_H,
                               facecolor='white', edgecolor='none', alpha=0.7, zorder=9))
        ax.text(span_left, label_y, iso_label, fontsize=8, ha='left', va='bottom', color=color, zorder=10)

    # 2) Reads
    order, drawn = (sorted(range(len(reads_blocks)), key=lambda i: reads_blocks[i][0][0]), set())
    rects, gap_rects, arrows = [], [], []
    for i in order:
        blks = reads_blocks[i]; intr = reads_introns[i]
        strand = read_strands[i]; rn = read_names[i]; y = disp_row(row_assign[i])
        grp = mapping.get(rn, '__UNASSIGNED__'); col = iso_colors.get(grp, ua_color)
        # Draw a thin continuous baseline from first to last block, then overlay blocks
        s_read, e_read = blks[0][0], blks[-1][1]
        gap_rects.append(Rectangle((s_read, y - INTRON_H/2), e_read - s_read, INTRON_H,
                                   facecolor=col, edgecolor='none'))
        for s, e in blks:
            rects.append(Rectangle((s, y - rh/2), e - s, rh, facecolor=col, edgecolor='none'))
        arr = '>' if strand == '+' else '<'
        arrows += [(blks[0][0], y, 'left', arr), (blks[-1][1], y, 'right', arr)]
    if rects:
        ax.add_collection(PatchCollection(rects, match_original=True))
    if gap_rects:
        ax.add_collection(PatchCollection(gap_rects, match_original=True))
    for x, y, ha, arr in arrows:
        ax.text(x, y, arr, ha=ha, va='center', fontsize=2, color='white')

    # 3) KDE/hist panels
    def _kde_or_hist(vals, lo, hi, n_bins=200):
        """Return left bin edges (xs) and normalized heights (ys in [0,1]).

        - Uses regularized, fixed-width bins across [lo, hi). This prevents
          the last bar from extending beyond the panel bounds.
        - For KDE, evaluate at bin centers and map to left edges for drawing.
        - Always normalize to max of 1 to limit vertical height within panel.
        """
        if not vals or hi <= lo:
            return None, None
        # Define edges and centers explicitly to avoid endpoint overshoot
        edges = np.linspace(lo, hi, n_bins + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        xs = edges[:-1]  # left edges for drawing
        try:
            from scipy.stats import gaussian_kde  # type: ignore
            if len(set(vals)) > 1:
                kde = gaussian_kde(vals)
                ys = kde(centers)
            else:
                ys = np.zeros_like(centers)
                # Put a single spike in the bin containing the value
                b = np.digitize([vals[0]], edges, right=False)[0] - 1
                b = max(0, min(b, len(ys) - 1))
                ys[b] = 1.0
        except Exception:
            counts, edges_h = np.histogram(vals, bins=n_bins, range=(lo, hi))
            ys = counts.astype(float)
            # Light smoothing for aesthetics
            if ys.size >= 5:
                kernel = np.ones(5) / 5.0
                ys = np.convolve(ys, kernel, mode='same')
        # Normalize safely
        mx = float(np.max(ys)) if np.size(ys) else 0.0
        if mx > 0:
            ys = ys / mx
        return xs, ys

    for iso in isos:
        if iso not in collapsed:
            continue
        label_r, br, reads_guard, tss_r, tts_r, spacer_r = collapsed[iso]
        if iso in ci_map:
            cs, ce = ci_map[iso]
        else:
            idxs = iso_to_idxs[iso]
            cs = min(reads_blocks[j][0][0] for j in idxs)
            ce = max(reads_blocks[j][-1][1] for j in idxs)
        if r0 is not None:
            cs = max(cs, r0); ce = min(ce, r1)
            if ce <= cs:
                continue
        blks_iso = iso_blocks.get(iso, [(cs, ce)])
        if r0 is not None:
            blks_iso = [(max(bs, r0), min(be, r1)) for (bs, be) in blks_iso if be > r0 and bs < r1]
            if not blks_iso:
                continue
        color = iso_colors.get(iso, 'black')
        idxs = iso_to_idxs.get(iso, [])
        tis, tts = [], []
        for j in idxs:
            rb = reads_blocks[j]
            if not rb:
                continue
            s_read, e_read = rb[0][0], rb[-1][1]
            if read_strands[j] == '+':
                tis.append(s_read); tts.append(e_read)
            else:
                tis.append(e_read); tts.append(s_read)
        span_left = min(bs for bs, _ in blks_iso); span_right = max(be for _, be in blks_iso)
        kde_h = KDE_H
        # Flip reserved rows for display
        tss_y = disp_row(tss_r)
        tts_y = disp_row(tts_r)
        xs, ys = _kde_or_hist(tis, span_left, span_right)
        if xs is not None and ys is not None and len(xs) > 0:
            bw = (xs[1] - xs[0]) if len(xs) > 1 else max(1.0, (span_right - span_left) / 200.0)
            # Panel background (keep a handle for clipping)
            tss_panel_rect = Rectangle((span_left, tss_y - kde_h/2), span_right - span_left, kde_h,
                                       facecolor='none', edgecolor='black', linewidth=0.5, zorder=3)
            ax.add_patch(tss_panel_rect)
            for xleft, yval in zip(xs, ys):
                if yval <= 0:
                    continue
                # Clamp horizontal extent to the panel
                left = max(span_left, xleft)
                width = min(bw, span_right - left)
                if width <= 0:
                    continue
                r = Rectangle((left, tss_y - kde_h/2), width, float(yval) * kde_h,
                              facecolor=color, edgecolor='none', alpha=0.45, zorder=4)
                r.set_clip_path(tss_panel_rect)
                ax.add_patch(r)
        xs2, ys2 = _kde_or_hist(tts, span_left, span_right)
        if xs2 is not None and ys2 is not None and len(xs2) > 0:
            bw2 = (xs2[1] - xs2[0]) if len(xs2) > 1 else max(1.0, (span_right - span_left) / 200.0)
            tts_panel_rect = Rectangle((span_left, tts_y - kde_h/2), span_right - span_left, kde_h,
                                       facecolor='none', edgecolor='black', linewidth=0.5, zorder=3)
            ax.add_patch(tts_panel_rect)
            for xleft, yval in zip(xs2, ys2):
                if yval <= 0:
                    continue
                left = max(span_left, xleft)
                width = min(bw2, span_right - left)
                if width <= 0:
                    continue
                r2 = Rectangle((left, tts_y - kde_h/2), width, float(yval) * kde_h,
                               facecolor=color, edgecolor='none', alpha=0.45, zorder=4)
                r2.set_clip_path(tts_panel_rect)
                ax.add_patch(r2)

    # No labels for unassigned reads

    # ---- Axes, limits, and polyA highlighting ----
    ax.set_yticks([])
    contig = chrom if chrom else (contig_for_title or "unknown")
    if r0 is not None:
        xmin, xmax = max(r0, read_min), min(r1, read_max)
        if xmax <= xmin:
            xmin, xmax = r0, r1
    else:
        xmin, xmax = read_min, read_max
    xpad = max((xmax - xmin) * 0.05, 50)
    ax.set_xlim(xmin - xpad, xmax + xpad)

    # Compute y-limits in display coordinates after flipping
    all_disp_rows = []
    all_disp_rows.extend([disp_row(r) for r in row_assign if r is not None])
    all_disp_rows.extend([disp_row(r) for rows in collapsed.values() for r in rows])
    y_lo = (min(all_disp_rows) if all_disp_rows else 0) - 2
    y_hi = (max(all_disp_rows) if all_disp_rows else 10) + 6

    if seq:
        window_size = 18
        min_as = 12
        left_x = r0 if r0 is not None else read_min
        base_width = 1
        regions = []
        for i in range(0, len(seq) - window_size + 1):
            window = seq[i:i + window_size]
            if window.count("A") >= min_as:
                start = left_x + i * base_width
                end = start + window_size * base_width
                regions.append((start, end))
        # Merge overlapping or adjacent windows to reduce patch count
        if regions:
            regions.sort(key=lambda x: x[0])
            merged = []
            cs, ce = regions[0]
            for s, e in regions[1:]:
                if s <= ce:
                    ce = max(ce, e)
                else:
                    merged.append((cs, ce))
                    cs, ce = s, e
            merged.append((cs, ce))
            # Draw vertical bands spanning the full plot height
            for start, end in merged:
                ax.axvspan(start, end, ymin=0, ymax=1, facecolor="red", alpha=0.12, zorder=0.5)
    else:
        logging.warning("[browser] No sequence available for polyA highlighting")

    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel(f"{contig}")
    ax.ticklabel_format(style='sci', axis='x', scilimits=(6,6))
    fig.subplots_adjust(bottom=0.18)

    # ---- Save ----
    if skip_plot:
        plt.close(fig)
        return None
    out_png = outdir / f"{genome}_{contig}_{xmin}-{xmax}.png"
    fig.savefig(out_png, dpi=fig.dpi)
    plt.close(fig)
    width_px, height_px = int(fig_w * fig.dpi), int(fig_h * fig.dpi)
    logger = logging.getLogger(__name__)
    logger.debug(f"Figure dimensions: {width_px} x {height_px} pixels")
    logger.debug(f"Saved: {out_png}")
    return out_png


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    cfg = Config.from_json(args.config)
    generate(cfg, args.region)


if __name__ == '__main__':
    main()
