#!/usr/bin/env python3
"""
flair_automate/slicer.py

Extract region-specific files for FLAIR pipeline, pulling inputs
from your FLAIR‐correct outputs:

  <outdir>/correct/<sample>/<corr_run>/<sample>.corrected.bed
  <outdir>/correct/<sample>/<corr_run>/<sample>.corrected.bam

Library API:
  slice_region(...)
  slice_all_regions(...)

CLI:
  python slicer.py --gtf … --junctions … --reads_fasta … --chrom … etc.
"""
import os
import shlex
import subprocess
import argparse
import glob
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

def run(cmd, stdout=None, cwd=None):
    """Run subprocess, merging stdout/stderr."""
    subprocess.run(cmd, stdout=stdout, stderr=subprocess.STDOUT, check=True, cwd=cwd)

def slice_bam(full_bam, chrom, start, end, out_bam, dry=False):
    if dry:
        logging.info("[DRY] samtools view -b %s %s:%d-%d > %s",
                     full_bam, chrom, start, end, out_bam)
        return
    if os.path.exists(out_bam):
        logging.info("  BAM slice exists (%s), skipping", out_bam)
    else:
        logging.info("  Slicing BAM %s:%d-%d → %s", chrom, start, end, out_bam)
        run(["samtools", "view", "-b", full_bam, f"{chrom}:{start}-{end}"],
            stdout=open(out_bam, "wb"))
        run(["samtools", "index", out_bam])

def awk_slice(in_path, chrom, start, end, out_path,
              col_start=2, col_end=3, dry=False):
    if dry:
        logging.info("[DRY] awk slice %s:%d-%d → %s",
                     in_path, chrom, start, end, out_path)
        return
    if os.path.exists(out_path):
        logging.info("  Slice exists (%s), skipping", out_path)
        return
    logging.info("  Slicing %s on %s:%d-%d → %s",
                 os.path.basename(in_path), chrom, start, end, out_path)
    awk_cmd = (
        f"awk '$1==\"{chrom}\" && ${col_start}>={start} && ${col_end}<={end}' "
        f"{shlex.quote(in_path)} > {shlex.quote(out_path)}"
    )
    run(["bash", "-c", awk_cmd])

def slice_gtf_fully_contained(gtf_path, chrom, start, end, out_path, dry=False):
    if dry:
        logging.info("[DRY] slice GTF %s:%d-%d → %s",
                     chrom, start, end, out_path)
        return
    if os.path.exists(out_path):
        logging.info("  GTF slice exists (%s), skipping", out_path)
        return
    logging.info("  Slicing GTF %s on %s:%d-%d → %s",
                 os.path.basename(gtf_path), chrom, start, end, out_path)
    import re
    target = chrom if chrom.startswith("chr") else f"chr{chrom}"
    feats, lines, genes = {}, {}, []
    with open(gtf_path) as fh:
        for L in fh:
            if L.startswith("#") or not L.strip(): continue
            cols = L.rstrip("\n").split("\t")
            if cols[0] != target: continue
            typ = cols[2]; st, en = int(cols[3]), int(cols[4])
            if typ == "gene":
                genes.append((st, en, L))
            m = re.search(r'transcript_id\s+"([^"]+)"', cols[8])
            if not m: continue
            tid = m.group(1)
            feats.setdefault(tid, []).append((st, en))
            lines.setdefault(tid, []).append(L)
    overlaps = [
        iv
        for ivs in feats.values()
        for iv in ivs
        if not (iv[1] < start or iv[0] > end)
    ]
    ns, ne = ((min(iv[0] for iv in overlaps), max(iv[1] for iv in overlaps))
              if overlaps else (start, end))
    with open(out_path, "w") as out:
        for g_st, g_en, gL in sorted(genes):
            if g_st >= ns and g_en <= ne:
                out.write(gL)
        for tid, ivs in feats.items():
            if all(iv[0] >= ns and iv[1] <= ne for iv in ivs):
                for L in sorted(lines[tid], key=lambda L: int(L.split("\t")[3])):
                    out.write(L)

def mk_region_fasta(bed_path, fa_full, fa_out, dry=False):
    if dry:
        logging.info("[DRY] mk_region_fasta %s from %s → %s",
                     bed_path, fa_full, fa_out)
        return
    if os.path.exists(fa_out):
        logging.info("  FASTA exists (%s), skipping", fa_out)
        return
    logging.info("  Building FASTA slice %s → %s",
                 os.path.basename(fa_full), os.path.basename(fa_out))
    ids = set()
    with open(bed_path) as fh:
        for L in fh:
            parts = L.rstrip("\n").split("\t")
            if len(parts) >= 4:
                ids.add(parts[3])
    write = False
    with open(fa_full) as inf, open(fa_out, "w") as outf:
        for L in inf:
            if L.startswith(">"):
                hdr = L[1:].split()[0]
                write = hdr in ids
            if write:
                outf.write(L)

def slice_region(full_bam, query_bed, gtf, junctions,
                 reads, exp5_bed, exp3_bed,
                 ref5_bed, ref3_bed,
                 chrom, start, end, outdir, dry=False):
    """
    Slice all inputs for one region into `outdir`.
    """
    os.makedirs(outdir, exist_ok=True)
    logging.info("→ slice_region into %s", outdir)

    # 1) BAM
    slice_bam(full_bam, chrom, start, end,
              os.path.join(outdir, f"{chrom}_{start}-{end}.bam"), dry)

    # 2) query BED
    awk_slice(query_bed, chrom, start, end,
              os.path.join(outdir, f"{chrom}_{start}-{end}.bed"), 2, 3, dry)

    # 3) experiment 5′/3′
    if exp5_bed:
        awk_slice(exp5_bed, chrom, start, end,
                  os.path.join(outdir, f"{chrom}_{start}-{end}.exp5.bed"), 2, 3, dry)
    if exp3_bed:
        awk_slice(exp3_bed, chrom, start, end,
                  os.path.join(outdir, f"{chrom}_{start}-{end}.exp3.bed"), 2, 3, dry)

    # 4) reference 5′/3′
    if ref5_bed:
        awk_slice(ref5_bed, chrom, start, end,
                  os.path.join(outdir, f"{chrom}_{start}-{end}.ref5.bed"), 2, 3, dry)
    if ref3_bed:
        awk_slice(ref3_bed, chrom, start, end,
                  os.path.join(outdir, f"{chrom}_{start}-{end}.ref3.bed"), 2, 3, dry)

    # 5) GTF
    slice_gtf_fully_contained(gtf, chrom, start, end,
              os.path.join(outdir, f"{chrom}_{start}-{end}.gtf"), dry)

    # 6) junctions
    if junctions:
        awk_slice(junctions, chrom, start, end,
                  os.path.join(outdir, f"{chrom}_{start}-{end}.tab"), 2, 3, dry)

    # 7) FASTA reads
    for fa in reads:
        mk_region_fasta(
            os.path.join(outdir, f"{chrom}_{start}-{end}.bed"),
            fa,
            os.path.join(outdir, f"{chrom}_{start}-{end}.fasta"),
            dry
        )

def slice_all_regions(cfg, outdir, dry_run=False):
    """
    For each region in cfg["regions"], and for each
    <sample>/<corr_run> under <outdir>/correct/, slice all needed files.
    """
    regions      = cfg["regions"]
    gtf          = cfg["gtf"]
    junctions    = cfg.get("junctions")
    reads        = cfg["reads_fasta"]
    exp5_bed     = cfg.get("experiment_5_prime_regions_bed_file")
    exp3_bed     = cfg.get("experiment_3_prime_regions_bed_file")
    ref5_bed     = cfg.get("reference_5_prime_regions_bed_file")
    ref3_bed     = cfg.get("reference_3_prime_regions_bed_file")

    correct_base = os.path.join(outdir, "correct")

    for chrom, spans in regions.items():
        for span in spans:
            start, end  = map(int, span.split(":",1))
            region_id   = f"{chrom}_{start}_{end}"
            region_root = os.path.join(outdir, "regions", region_id)

            # find every corrected run
            pattern = os.path.join(correct_base, "*", "*")
            for sample_run in sorted(glob.glob(pattern)):
                sample   = os.path.basename(os.path.dirname(sample_run))
                corr_run = os.path.basename(sample_run)

                # locate primary outputs
                bam_files = glob.glob(os.path.join(sample_run, f"{sample}*.bam"))
                bed_files = glob.glob(os.path.join(sample_run, f"{sample}*.bed"))
                if not bam_files or not bed_files:
                    logging.warning("No BAM/BED in %s, skipping", sample_run)
                    continue

                full_bam  = bam_files[0]
                query_bed = bed_files[0]
                outdir_r  = os.path.join(region_root, sample, corr_run)

                slice_region(
                    full_bam, query_bed, gtf, junctions,
                    reads, exp5_bed, exp3_bed,
                    ref5_bed, ref3_bed,
                    chrom, start, end, outdir_r, dry_run
                )

# ─── CLI for ad-hoc slicing ────────────────────────────────────────────────────
if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Slice one region for FLAIR")
    p.add_argument('--gtf',       required=True)
    p.add_argument('--junctions', help="optional")
    p.add_argument('--reads_fasta', action='append', required=True)
    p.add_argument('--experiment_5_prime_regions_bed_file', help="optional")
    p.add_argument('--experiment_3_prime_regions_bed_file', help="optional")
    p.add_argument('--reference_5_prime_regions_bed_file',  help="optional")
    p.add_argument('--reference_3_prime_regions_bed_file',  help="optional")
    p.add_argument('--chrom',     required=True)
    p.add_argument('--start',     type=int, required=True)
    p.add_argument('--end',       type=int, required=True)
    p.add_argument('--outdir',    required=True)
    p.add_argument('--dry_run',   action='store_true')
    args = p.parse_args()

    slice_region(
        full_bam    = args.bam,
        query_bed   = args.query,
        gtf         = args.gtf,
        junctions   = args.junctions,
        reads       = args.reads_fasta,
        exp5_bed    = args.experiment_5_prime_regions_bed_file,
        exp3_bed    = args.experiment_3_prime_regions_bed_file,
        ref5_bed    = args.reference_5_prime_regions_bed_file,
        ref3_bed    = args.reference_3_prime_regions_bed_file,
        chrom       = args.chrom,
        start       = args.start,
        end         = args.end,
        outdir      = args.outdir,
        dry         = args.dry_run
    )




