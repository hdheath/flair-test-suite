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

import os
import json
import argparse
import warnings
from collections import defaultdict
from itertools import chain

import pandas as pd
import pysam
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from scipy.stats import gaussian_kde
from intervaltree import IntervalTree
import time
from matplotlib.collections import PatchCollection, LineCollection


def parse_args():
    p = argparse.ArgumentParser(
        description="Reads/genes/collapsed isoforms + reversed read-span above"
    )
    p.add_argument('-c','--config', required=True,
        help="JSON config keys: bam, gtf, mapping, collapsed_isoforms, genome, outdir, "
             "gene_height, gene_row_height, read_row_height, fig_width")
    return p.parse_args()


def load_gtf(gtf_path):
    if not os.path.exists(gtf_path):
        warnings.warn(f"GTF not found: {gtf_path}")
        return pd.DataFrame()
    df = pd.read_csv(
        gtf_path, sep='\t', header=None, comment='#',
        names=['chr','source','feature','start','end','score','strand','frame','attributes']
    )
    tx = df[df['feature']=='transcript'].copy()
    tx['gene_id'] = tx['attributes'].str.extract(r'gene_id "([^"]+)"')
    collapsed = (
        tx
        .groupby(['gene_id','strand'])
        .agg({'chr':'first','start':'min','end':'max'})
        .reset_index()
    )
    return collapsed[['chr','start','end','strand','gene_id']]


def load_mapping(map_path):
    m = {}
    if not map_path or not os.path.exists(map_path):
        warnings.warn(f"Mapping not found: {map_path}")
        return m
    with open(map_path) as f:
        for line in f:
            iso, reads = line.strip().split('\t')
            for r in reads.split(','):
                if r:
                    m[r] = iso
    if not m:
        warnings.warn("Mapping is empty.")
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


def main():
    args = parse_args()
    cfg = json.load(open(args.config))

    # timing start
    t_global_start = time.time()

    # config
    bam          = cfg['bam']
    gtf          = cfg['gtf']
    mapping_file = cfg.get('mapping')
    cis_bed      = cfg.get('collapsed_isoforms')
    genome       = cfg.get('genome', '')
    outdir       = cfg.get('outdir', '.')
    gene_h       = cfg.get('gene_height', 1.5)
    gene_row_h   = cfg.get('gene_row_height', 3)
    read_row_h   = cfg.get('read_row_height', 0.02)
    fig_w        = cfg.get('fig_width', 12)
    os.makedirs(outdir, exist_ok=True)

    # load reads
    reads_blocks, reads_introns, read_strands, read_names = [], [], [], []
    bf = pysam.AlignmentFile(bam, 'rb')
    for rd in bf.fetch():
        if rd.is_unmapped: continue
        blks = rd.get_blocks()
        if not blks: continue
        reads_blocks.append(blks)
        reads_introns.append([(blks[i][1], blks[i+1][0]) for i in range(len(blks)-1)])
        read_strands.append('-' if rd.is_reverse else '+')
        read_names.append(rd.query_name)
    if not reads_blocks:
        warnings.warn("No reads loaded.")
        return

    all_starts = [b[0] for blks in reads_blocks for b in blks]
    all_ends   = [b[1] for blks in reads_blocks for b in blks]
    read_min, read_max = min(all_starts), max(all_ends)

    # load mapping & collapsed isoforms
    mapping = load_mapping(mapping_file)
    cis_df = pd.DataFrame(columns=['iso_id','start','end','strand'])
    if cis_bed and os.path.exists(cis_bed):
        tmp = pd.read_csv(cis_bed, sep='\t', header=None, comment='#').iloc[:,:6]
        tmp.columns = ['chr','start','end','iso_id','score','strand']
        cis_df = tmp[['iso_id','start','end','strand']]
    ci_map = {r.iso_id:(r.start, r.end) for _,r in cis_df.iterrows()}

    # group by isoform
    iso_to_idxs, unassigned = defaultdict(list), []
    for i,name in enumerate(read_names):
        iso = mapping.get(name)
        if iso: iso_to_idxs[iso].append(i)
        else:   unassigned.append(i)

    palette = list(chain(plt.cm.tab20.colors, plt.cm.Set3.colors, plt.cm.Dark2.colors))
    pal = [c for c in palette if max(c)<0.9 and not (max(c)-min(c)<0.2 and max(c)>0.7)]
    isos = sorted(iso_to_idxs, key=lambda i: len(iso_to_idxs[i]), reverse=True)
    iso_colors = {iso: pal[i%len(pal)] for i,iso in enumerate(isos)}
    ua_color = 'lightgrey'

    # layout
    occupancy, row_assign, collapsed = [], [None]*len(reads_blocks), {}
    ih, rh, new_h = gene_h*0.6, gene_h*0.6, gene_h*3

    for iso in isos:
        idxs = iso_to_idxs[iso]
        sorted_idxs = sorted(idxs, key=lambda i: reads_blocks[i][0][0])
        blks_list,_ = zip(*[(reads_blocks[i],None) for i in sorted_idxs])
        rel,_ = assign_read_rows(blks_list)

        st,en = ci_map.get(iso,(None,None))
        if st is None:
            warnings.warn(f"No collapsed bar for {iso}")
            continue

        left = min(reads_blocks[j][0][0] for j in idxs)
        right= max(reads_blocks[j][-1][1] for j in idxs)
        max_rel = max(rel)
        specs = [(1,st,en),(2,st,en),(4,left,right),(6,left,right),(9,left,right),(11,left,right),(14,left,right)]

        for base in range(len(occupancy)+1):
            ok=True
            for ri,blks in zip(rel,blks_list):
                for off in (0,1):
                    tgt=base+ri+off
                    if tgt<len(occupancy) and any(len(occupancy[tgt].overlap(b[0],b[1]))>0 for b in blks): ok=False;break
                if not ok: break
            if not ok: continue
            for off,ps,pe in specs:
                pr=base+max_rel+off
                if pr<len(occupancy) and occupancy[pr].overlap(ps,pe): ok=False; break
            if ok: break
        else:
            base=len(occupancy)

        max_r=base+max_rel; space1=max_r+1; col_r=max_r+2; space2=max_r+4
        emp1 =max_r+6; space3=max_r+9; emp2=max_r+11; space4=max_r+14
        for r in range(space4+1):
            if r>=len(occupancy): occupancy.append(IntervalTree())
        for ri,idx in zip(rel,sorted_idxs):
            r=base+ri; row_assign[idx]=r
            for b in reads_blocks[idx]: occupancy[r].addi(b[0],b[1],True)
        for r in (space1,col_r,space2,emp1,space3,emp2,space4):
            occupancy[r].addi(left if r>=space2 else st, right if r>=space2 else en, True)
        collapsed[iso]=(space1,col_r,space2,emp1,space3,emp2,space4)

    for idx in unassigned:
        blks=reads_blocks[idx]; placed=False
        for r in range(len(occupancy)):
            if any(len(occupancy[r].overlap(b[0],b[1]))>0 for b in blks): continue
            if r>0 and any(len(occupancy[r-1].overlap(b[0],b[1]))>0 for b in blks): continue
            if r<len(occupancy)-1 and any(len(occupancy[r+1].overlap(b[0],b[1]))>0 for b in blks): continue
            row_assign[idx]=r
            for b in blks: occupancy[r].addi(b[0],b[1],True)
            placed=True; break
        if not placed:
            t=IntervalTree()
            for b in blks: t[b[0]:b[1]]=True
            occupancy.append(t); row_assign[idx]=len(occupancy)-1
            warnings.warn(f"Unassigned read at row {row_assign[idx]}")

    tot=len(occupancy)
    row_assign=[tot-1-r for r in row_assign]
    collapsed={iso:tuple(tot-1-r for r in rows) for iso,rows in collapsed.items()}

    gtf_df=load_gtf(gtf)
    if gtf_df.empty:
        warnings.warn("No genes.");return
    contig,rs,re_=(gtf_df.chr.iloc[0],int(gtf_df.start.min()),int(gtf_df.end.max()))
    gb_buffered=[(max(0,s-5),e+5) for s,e in zip(gtf_df.start,gtf_df.end)]
    rg,gn=assign_read_rows([[b] for b in gb_buffered])
    lab, gsp = gene_h*0.3, gene_h*7
    gbase=max(row_assign)+2
    gene_y={i:gbase+rg[i]*gsp for i in range(len(rg))}

    fig_h=read_row_h*tot+gene_row_h*gn
    fig=plt.figure(figsize=(fig_w,fig_h))
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
            # KDE
            tis,tts=[],[]
            for j in idxs:
                s,e=reads_blocks[j][0][0],reads_blocks[j][-1][1]
                if read_strands[j]=='+': tis.append(s); tts.append(e)
                else: tis.append(e); tts.append(s)
            bins=np.linspace(left,right,300)
            bw=bins[1]-bins[0]
            if len(set(tis))>1:
                h=gaussian_kde(tis)(bins); h/=h.max()
                for xval,hval in zip(bins,h): ax.add_patch(Rectangle((xval,er-new_h/2),bw,hval*new_h,facecolor=iso_colors[iso],edgecolor='none',alpha=0.6,zorder=4))
            if len(set(tts))>1:
                h=gaussian_kde(tts)(bins); h/=h.max()
                for xval,hval in zip(bins,h): ax.add_patch(Rectangle((xval,e2-new_h/2),bw,hval*new_h,facecolor=iso_colors[iso],edgecolor='none',alpha=0.6,zorder=4))
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
    smoothed_tts = np.convolve(counts_tts, np.ones(5)/5, mode='same')
    if smoothed_tts.max()>0:
        smoothed_tts /= smoothed_tts.max()
    width_bin = edges_tts[1] - edges_tts[0]
    for left_edge, h in zip(edges_tts[:-1], smoothed_tts):
        ax.add_patch(Rectangle(
            (left_edge, y1), width_bin, h * rect_h,
            facecolor='red', edgecolor='none', alpha=0.4, zorder=3
        ))

    # top box: plot TSS histogram
    counts_tss, edges_tss = np.histogram(tss_pos, bins=200, weights=tss_w)
    smoothed_tss = np.convolve(counts_tss, np.ones(5)/5, mode='same')
    if smoothed_tss.max()>0:
        smoothed_tss /= smoothed_tss.max()
    width_bin = edges_tss[1] - edges_tss[0]
    for left_edge, h in zip(edges_tss[:-1], smoothed_tss):
        ax.add_patch(Rectangle(
            (left_edge, y2), width_bin, h * rect_h,
            facecolor='blue', edgecolor='none', alpha=0.4, zorder=3
        ))


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
    out_png = os.path.join(outdir, f"{genome}_{contig}_{rs}-{re_}.png")
    plt.savefig(out_png, dpi=600, bbox_inches='tight')
    plt.close(fig)
    # print figure dimensions in pixels
    width_in, height_in = fig.get_size_inches()
    dpi = 600
    width_px, height_px = int(width_in * dpi), int(height_in * dpi)
    print(f"Figure dimensions: {width_px} x {height_px} pixels")
    print(f"Saved: {out_png}")


if __name__ == '__main__':
    main()