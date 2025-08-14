#!/usr/bin/env python3
"""
plot_genome_browser_v13.py

Generate a genome-browser style plot showing gene annotations and aligned reads with splice junctions.
Takes a JSON config file as input with keys:
  bam: path to BAM
  gtf: path to GTF
  junctions: path to splice junctions .tab (optional)
  mapping: path to isoform–read mapping file (optional)
  genome: genome name (for title)
  outdir: output directory
  gene_height: vertical size of gene exon rectangles [default: 1.5]
  gene_row_height: vertical size per gene row [default: 3]
  read_row_height: vertical size per read row [default: 0.02]
  fig_width: figure width in inches [default: 12]

Features:
1) Isoforms sorted by descending read count—each isoform’s reads plotted as a block before moving to next.
2) Unassigned reads placed individually: for each, find the last isoform row overlapping its x-range, then place immediately below that (or next available row).
3) Group-first-fit for isoform blocks.
4) Avoid near-white/grey isoform colors by filtering colormaps.
5) After assignment, rows flipped so largest-coverage isoform at top and unassigned trailing.
6) Tight vertical spacing (hspace=0) and fixed margins above/below reads.

python ted_browser_v2.4.py -c /private/groups/brookslab/hdheath/projects/ted/transcriptome_browser/ted_browser_test_v.3.json

"""

import json
import argparse
import logging
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
from scipy.stats import gaussian_kde
from intervaltree import IntervalTree
from matplotlib.collections import PatchCollection, LineCollection


@dataclass(frozen=True)
class Config:
    bam: Path
    gtf: Path
    genome: str = ""
    outdir: Path = Path(".")
    mapping: Optional[Path] = None
    collapsed_isoforms: Optional[Path] = None
    gene_height: float = 1.5
    gene_row_height: float = 3.0
    read_row_height: float = 0.02
    fig_width: float = 12.0
    iso_kde: bool = True

    @staticmethod
    def from_json(p: Path) -> "Config":
        cfg = json.loads(p.read_text())
        return Config(
            bam=Path(cfg["bam"]),
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
        )


def parse_args():
    p = argparse.ArgumentParser(
        description="Reads/genes/collapsed isoforms + reversed read-span above",
    )
    p.add_argument(
        "-c",
        "--config",
        type=Path,
        required=True,
        help=(
            "JSON config keys: bam, gtf, mapping, collapsed_isoforms, genome, outdir, "
            "gene_height, gene_row_height, read_row_height, fig_width, iso_kde"
        ),
    )
    p.add_argument(
        "--region",
        help="Limit plotting region (e.g. chr1:100-200)",
    )
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


def reserve_iso_rows(base: int, max_rel: int) -> dict[str, int]:
    return {
        "space1": base + max_rel + 1,
        "collapsed_bar": base + max_rel + 2,
        "space2": base + max_rel + 4,
        "empty1": base + max_rel + 6,
        "space3": base + max_rel + 9,
        "empty2": base + max_rel + 11,
        "space4": base + max_rel + 14,
    }


def rects_from_bins(edges, heights, y0, hscale, facecolor):
    rects = []
    width = edges[1] - edges[0]
    for x, h in zip(edges[:-1], heights):
        if h <= 0:
            continue
        rects.append(Rectangle((x, y0), width, h * hscale))
    pc = PatchCollection(rects, match_original=False)
    pc.set_facecolor(facecolor)
    pc.set_edgecolor("none")
    pc.set_alpha(0.4)
    pc.set_zorder(3)
    return pc


def load_gtf(gtf_path: Path, region=None):
    if not gtf_path.exists():
        logging.warning(f"GTF not found: {gtf_path}")
        return pd.DataFrame()
    df = pd.read_csv(
        gtf_path,
        sep='\t',
        header=None,
        comment='#',
        names=['chr','source','feature','start','end','score','strand','frame','attributes'],
    )
    tx = df[df['feature']=='transcript'].copy()
    tx['gene_id'] = tx['attributes'].str.extract(r'gene_id "([^"]+)"')
    collapsed = (
        tx.groupby(['gene_id','strand']).agg({'chr':'first','start':'min','end':'max'}).reset_index()
    )
    res = collapsed[['chr','start','end','strand','gene_id']]
    if region:
        chrom, start, end = region
        res = res[(res['chr']==chrom) & (res['start']<end) & (res['end']>start)]
    return res


def load_mapping(map_path: Optional[Path]):
    m: dict[str, str] = {}
    if not map_path or not map_path.exists():
        logging.warning(f"Mapping not found: {map_path}")
        return m
    with map_path.open() as f:
        for line in f:
            iso, reads = line.strip().split('\t')
            for r in reads.split(','):
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

def _parse_region(region: Optional[str]) -> Optional[Tuple[str, int, int]]:
    """Parse a region string like 'chr1:100-200' into a tuple."""
    if not region:
        return None
    try:
        chrom, start, end = region.replace(":", "\t").replace("-", "\t").split("\t")
        start_i, end_i = int(start), int(end)
    except Exception:
        logging.warning(f"Could not parse region string: {region}")
        return None
    if end_i - start_i >= 20000:
        logging.warning("Region length >= 20000bp; skipping plot")
        return None
    return chrom, start_i, end_i


def generate(cfg: Config, region: Optional[str] = None) -> None:
    """Generate transcriptome browser plot from configuration."""
    region_tuple = _parse_region(region)

    bam = cfg.bam
    gtf = cfg.gtf
    if not bam.exists():
        logging.warning(f"BAM not found: {bam}")
        return
    if not gtf.exists():
        logging.warning(f"GTF not found: {gtf}")
        return
    mapping_file = cfg.mapping if cfg.mapping and cfg.mapping.exists() else None
    if cfg.mapping and not mapping_file:
        logging.warning(f"Mapping not found: {cfg.mapping}")
    cis_bed = cfg.collapsed_isoforms if cfg.collapsed_isoforms and cfg.collapsed_isoforms.exists() else None
    if cfg.collapsed_isoforms and not cis_bed:
        logging.warning(f"Collapsed isoforms BED not found: {cfg.collapsed_isoforms}")
    genome = cfg.genome
    outdir = cfg.outdir
    gene_h = cfg.gene_height
    gene_row_h = cfg.gene_row_height
    read_row_h = cfg.read_row_height
    fig_w = cfg.fig_width
    outdir.mkdir(parents=True, exist_ok=True)

    # load reads
    reads_blocks, reads_introns, read_strands, read_names = [], [], [], []
    with pysam.AlignmentFile(bam, "rb") as bf:
        iterator = bf.fetch(*region_tuple) if region_tuple else bf.fetch()
        for rd in iterator:
            if rd.is_unmapped:
                continue
            blks = rd.get_blocks()
            if not blks:
                continue
            reads_blocks.append(blks)
            reads_introns.append([(blks[i][1], blks[i + 1][0]) for i in range(len(blks) - 1)])
            read_strands.append(read_strand(rd))
            read_names.append(rd.query_name)
    if not reads_blocks:
        logging.warning("No reads loaded.")
        return

    all_starts = [b[0] for blks in reads_blocks for b in blks]
    all_ends = [b[1] for blks in reads_blocks for b in blks]
    read_min, read_max = min(all_starts), max(all_ends)

    # load mapping & collapsed isoforms
    mapping = load_mapping(mapping_file)
    cis_df = pd.DataFrame(columns=['iso_id','start','end','strand'])
    if cis_bed and cis_bed.exists():
        tmp = pd.read_csv(cis_bed, sep='\t', header=None, comment='#').iloc[:, :6]
        tmp.columns = ['chr','start','end','iso_id','score','strand']
        cis_df = tmp[['iso_id','start','end','strand']]
    ci_map = {r.iso_id: (r.start, r.end) for _, r in cis_df.iterrows()}

    # group by isoform
    iso_to_idxs, unassigned = defaultdict(list), []
    for i,name in enumerate(read_names):
        iso = mapping.get(name)
        if iso: iso_to_idxs[iso].append(i)
        else:   unassigned.append(i)

    base_colors = list(chain(plt.cm.tab20.colors, plt.cm.Set3.colors, plt.cm.Dark2.colors))
    pal = pleasant_colors(base_colors)
    isos = sorted(iso_to_idxs, key=lambda i: len(iso_to_idxs[i]), reverse=True)
    iso_colors = {iso: pal[i%len(pal)] for i,iso in enumerate(isos)}
    ua_color = 'lightgrey'

    # layout
    occupancy, row_assign, collapsed = [], [None] * len(reads_blocks), {}
    ih, rh, new_h = gene_h * 0.6, gene_h * 0.6, gene_h * 3

    for iso in isos:
        idxs = iso_to_idxs[iso]
        sorted_idxs = sorted(idxs, key=lambda i: reads_blocks[i][0][0])
        blks_list, _ = zip(*[(reads_blocks[i], None) for i in sorted_idxs])
        rel, _ = assign_read_rows(blks_list)

        st, en = ci_map.get(iso, (None, None))
        if st is None:
            logging.warning(f"No collapsed bar for {iso}")
            continue

        left = min(reads_blocks[j][0][0] for j in idxs)
        right = max(reads_blocks[j][-1][1] for j in idxs)
        max_rel = max(rel)

        for base in range(len(occupancy) + 1):
            ok = True
            for ri, blks in zip(rel, blks_list):
                for off in (0, 1):
                    tgt = base + ri + off
                    if tgt < len(occupancy) and overlaps(occupancy[tgt], blks):
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                continue
            spec_rows = reserve_iso_rows(base, max_rel)
            for key, pr in spec_rows.items():
                seg = (left, right) if key not in ("space1", "collapsed_bar") else (st, en)
                if pr < len(occupancy) and overlaps(occupancy[pr], [seg]):
                    ok = False
                    break
            if ok:
                break
        else:
            base = len(occupancy)
            spec_rows = reserve_iso_rows(base, max_rel)

        ensure_rows(occupancy, spec_rows["space4"])
        for ri, idx in zip(rel, sorted_idxs):
            r = base + ri
            row_assign[idx] = r
            add_blocks(occupancy[r], reads_blocks[idx])
        for key in ("space1", "collapsed_bar", "space2", "empty1", "space3", "empty2", "space4"):
            r = spec_rows[key]
            seg = (left, right) if key not in ("space1", "collapsed_bar") else (st, en)
            add_blocks(occupancy[r], [seg])
        collapsed[iso] = (
            spec_rows["space1"],
            spec_rows["collapsed_bar"],
            spec_rows["space2"],
            spec_rows["empty1"],
            spec_rows["space3"],
            spec_rows["empty2"],
            spec_rows["space4"],
        )

    for idx in unassigned:
        blks = reads_blocks[idx]
        placed = False
        for r in range(len(occupancy)):
            if overlaps(occupancy[r], blks):
                continue
            if r > 0 and overlaps(occupancy[r - 1], blks):
                continue
            if r < len(occupancy) - 1 and overlaps(occupancy[r + 1], blks):
                continue
            row_assign[idx] = r
            add_blocks(occupancy[r], blks)
            placed = True
            break
        if not placed:
            t = IntervalTree()
            add_blocks(t, blks)
            occupancy.append(t)
            row_assign[idx] = len(occupancy) - 1
            logging.warning(f"Unassigned read at row {row_assign[idx]}")

    tot=len(occupancy)
    row_assign=[tot-1-r for r in row_assign]
    collapsed={iso:tuple(tot-1-r for r in rows) for iso,rows in collapsed.items()}

    gtf_df=load_gtf(gtf, region_tuple)
    if gtf_df.empty:
        logging.warning("No genes.")
        return
    contigs = gtf_df["chr"].unique().tolist()
    if len(contigs) != 1:
        logging.warning(f"Multiple contigs in GTF view: {contigs}. Using first for title.")
    contig = contigs[0]
    rs, re_ = int(gtf_df.start.min()), int(gtf_df.end.max())
    gb_buffered=[(max(0,s-5),e+5) for s,e in zip(gtf_df.start,gtf_df.end)]
    rg,gn=assign_read_rows([[b] for b in gb_buffered])
    lab, gsp = gene_h*0.3, gene_h*7
    gbase=max(row_assign)+2
    gene_y={i:gbase+rg[i]*gsp for i in range(len(rg))}

    fig_h=read_row_h*tot+gene_row_h*gn
    fig=plt.figure(figsize=(fig_w,fig_h), dpi=600)
    ax=fig.add_subplot(1,1,1)
    for side in ['left','right','bottom','top']:
        ax.spines[side].set_color('gray'); ax.spines[side].set_linewidth(0.4); ax.spines[side].set_alpha(0.4)

    order, drawn = (sorted(range(len(reads_blocks)), key=lambda i: reads_blocks[i][0][0]),set())
    rects, lines, arrows=[],[],[]
    for i in order:
        blks=reads_blocks[i]; intr=reads_introns[i]
        strand=read_strands[i]; rn=read_names[i]; y=row_assign[i]
        grp=mapping.get(rn,'__UNASSIGNED__'); col=iso_colors.get(grp,ua_color)
        for s,e in blks: rects.append(Rectangle((s,y-rh/2),e-s,rh,facecolor=col,linewidth=0))
        for a,b in intr: lines.append([(a,y),(b,y)])
        arr='>' if strand=='+' else '<'
        arrows += [(blks[0][0],y,'left',arr),(blks[-1][1],y,'right',arr)]
    ax.add_collection(PatchCollection(rects,match_original=True))
    ax.add_collection(LineCollection(lines,colors='black',linewidths=0.5))
    for x,y,ha,arr in arrows: ax.text(x,y,arr,ha=ha,va='center',fontsize=2,color='white')

    # collapsed placeholders & per-iso KDE
    for iso in isos:
        if iso not in drawn:
            s1,br,s2,er,s3,e2,s4 = collapsed[iso]
            cs,ce=ci_map[iso]
            w=ce-cs
            idxs=iso_to_idxs[iso]
            left=min(reads_blocks[j][0][0] for j in idxs)
            right=max(reads_blocks[j][-1][1] for j in idxs)
            ax.add_patch(Rectangle((cs,s1-rh/2),w,rh,facecolor='none',edgecolor='none',zorder=1))
            ax.add_patch(Rectangle((cs,br-ih/2),w,ih,facecolor=iso_colors[iso],edgecolor='black',linewidth=0.1,zorder=3))
            ax.add_patch(Rectangle((left,s2-new_h/2),right-left,new_h*5,facecolor='none',edgecolor='none',zorder=2))
            ax.add_patch(Rectangle((left,er-new_h/2),right-left,new_h,facecolor='none',edgecolor='black',linewidth=0.1,zorder=2))
            ax.add_patch(Rectangle((left,s3-new_h/2),right-left,new_h,facecolor='none',edgecolor='none',zorder=2))
            ax.add_patch(Rectangle((left,e2-new_h/2),right-left,new_h,facecolor='none',edgecolor='black',linewidth=0.1,zorder=2))
            ax.add_patch(Rectangle((left,s4-new_h/2),right-left,new_h,facecolor='none',edgecolor='none',zorder=2))
            if cfg.iso_kde and len(idxs) >= 5:
                tis, tts = [], []
                for j in idxs:
                    s, e = reads_blocks[j][0][0], reads_blocks[j][-1][1]
                    if read_strands[j] == '+':
                        tis.append(s); tts.append(e)
                    else:
                        tis.append(e); tts.append(s)
                bins = np.linspace(left, right, 300)
                bw = bins[1] - bins[0]
                if len(set(tis)) > 1:
                    h = gaussian_kde(tis)(bins)
                    h /= h.max()
                    for xval, hval in zip(bins, h):
                        ax.add_patch(Rectangle((xval, er - new_h/2), bw, hval * new_h,
                                               facecolor=iso_colors[iso], edgecolor='none', alpha=0.6, zorder=4))
                if len(set(tts)) > 1:
                    h = gaussian_kde(tts)(bins)
                    h /= h.max()
                    for xval, hval in zip(bins, h):
                        ax.add_patch(Rectangle((xval, e2 - new_h/2), bw, hval * new_h,
                                               facecolor=iso_colors[iso], edgecolor='none', alpha=0.6, zorder=4))
            label_x=left-(right-left)*0.005
            ax.text(label_x,er,'TSS',ha='right',va='center',rotation=90,fontsize=3,color='black',alpha=0.7,zorder=5)
            ax.text(label_x,e2,'TTS',ha='right',va='center',rotation=90,fontsize=3,color='black',alpha=0.7,zorder=5)
            drawn.add(iso)

    # draw genes
    for i,row in gtf_df.iterrows():
        s,e=row.start,row.end; y=gene_y[i]
        c='blue' if row.strand=='+' else 'red'
        ax.add_patch(Rectangle((s,y-gene_h/2),e-s,gene_h,color=c,linewidth=1))
        step=max(150,(e-s)//10)
        arr='>' if row.strand=='+' else '<'
        for x in range(s+step,e-step,step): ax.text(x,y,arr,ha='center',va='center',fontsize=6,color='white')
        ax.text((s+e)/2,y+gene_h/2+lab,row.gene_id,ha='center',va='bottom',fontsize=8)

    # labels & limits
    ax.set_yticks([])
    ax.set_xlabel(f"Position on {contig}")
    ax.set_title(f"{genome}:{contig}:{rs}-{re_}",fontsize=14,pad=6)

    # ── compute limits & weighted KDE ──
    all_rows = row_assign + [r for rows in collapsed.values() for r in rows]
    y_lo     = min(all_rows) - new_h/2 - 1
    old_y_hi = max(gene_y.values()) + gene_h/2 + lab + 10
    x_lo     = rs - (re_ - rs)*0.05
    x_hi     = re_ + (re_ - rs)*0.05

    # build per-isoform weighted TSS/TTS
    tss_pos, tss_w = [], []
    tts_pos, tts_w = [], []
    for iso, idxs in iso_to_idxs.items():
        n = len(idxs)
        if n==0: continue
        w = 1.0 / n
        for j in idxs:
            blks = reads_blocks[j]
            s, e = blks[0][0], blks[-1][1]
            strand = read_strands[j]
            if strand=='+':
                tss_pos.append(s); tss_w.append(w)
                tts_pos.append(e); tts_w.append(w)
            else:
                tss_pos.append(e); tss_w.append(w)
                tts_pos.append(s); tts_w.append(w)
    total_tss = sum(tss_w)
    total_tts = sum(tts_w)
    if total_tss>0: tss_w = [wi/total_tss for wi in tss_w]
    if total_tts>0: tts_w = [wi/total_tts for wi in tts_w]

    # empty-rectangle geometry
    rect_h          = 8 * gene_h
    spacing_to_gene = 10 * gene_h
    spacing_between = 1  * gene_h
    y_gene_top      = max(gene_y.values()) + gene_h/2
    y1 = y_gene_top + spacing_to_gene
    y2 = y1 + rect_h + spacing_between
    width = read_max - read_min

    # draw empty rectangles
    ax.add_patch(Rectangle((read_min, y1), width, rect_h,
                           facecolor='none', edgecolor='black', linewidth=0.5, zorder=2))
    ax.add_patch(Rectangle((read_min, y2), width, rect_h,
                           facecolor='none', edgecolor='black', linewidth=0.5, zorder=2))

    # ── binned histogram + smoothing in bottom box (TSS) ──
    # bottom box: plot TTS histogram
    counts_tts, edges_tts = np.histogram(tts_pos, bins=200, weights=tts_w)
    smoothed_tts = np.convolve(counts_tts, np.ones(5) / 5, mode='same')
    if smoothed_tts.max() > 0:
        smoothed_tts /= smoothed_tts.max()
    ax.add_collection(rects_from_bins(edges_tts, smoothed_tts, y1, rect_h, 'red'))

    # top box: plot TSS histogram
    counts_tss, edges_tss = np.histogram(tss_pos, bins=200, weights=tss_w)
    smoothed_tss = np.convolve(counts_tss, np.ones(5) / 5, mode='same')
    if smoothed_tss.max() > 0:
        smoothed_tss /= smoothed_tss.max()
    ax.add_collection(rects_from_bins(edges_tss, smoothed_tss, y2, rect_h, 'blue'))


    # external labels: TTS on bottom, TSS on top
    x_lbl = read_min - width * 0.001
    y1c  = y1 + rect_h/2   # center of bottom box
    y2c  = y2 + rect_h/2   # center of top box

    ax.text(
        x_lbl, y1c, 'TTS',
        ha='right', va='center', rotation=90,
        fontsize=7, color='red', alpha=0.7, zorder=5
    )
    ax.text(
        x_lbl, y2c, 'TSS',
        ha='right', va='center', rotation=90,
        fontsize=7, color='blue', alpha=0.7, zorder=5
    )

    # finalize y_hi & x limits
    rect_top = y2 + rect_h
    y_hi     = max(old_y_hi, rect_top + gene_h)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)


    # save
    out_png = outdir / f"{genome}_{contig}_{rs}-{re_}.png"
    fig.savefig(out_png, bbox_inches='tight')
    plt.close(fig)
    width_px, height_px = int(fig_w * fig.dpi), int(fig_h * fig.dpi)
    print(f"Figure dimensions: {width_px} x {height_px} pixels")
    print(f"Saved: {out_png}")


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    cfg = Config.from_json(args.config)
    generate(cfg, args.region)


if __name__ == '__main__':
    main()
