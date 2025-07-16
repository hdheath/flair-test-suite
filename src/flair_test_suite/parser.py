"""
src/flair_test_suite/parser.py

Region-based slicing for FLAIR pipeline: extract per-region BAM, BED, GTF,
junctions, and FASTA slices based on `cfg.whole_sample`, `cfg.run.regions`,
and `cfg.sample_name`, placing outputs under `outputs/<sample_name>/…`.
"""

import logging
import subprocess
import shlex
import glob
from pathlib import Path


def run(cmd, stdout=None, cwd=None):
    """Run a subprocess, merging stdout/stderr."""
    subprocess.run(cmd, stdout=stdout, stderr=subprocess.STDOUT, check=True, cwd=cwd)


def slice_bam(full_bam: Path, chrom: str, start: int, end: int,
              out_bam: Path, dry: bool = False):
    """Extract a BAM slice for a given region."""
    if dry:
        logging.info("[DRY] samtools view -b %s %s:%d-%d > %s",
                     full_bam, chrom, start, end, out_bam)
        return
    if out_bam.exists():
        logging.info("BAM slice exists (%s), skipping", out_bam)
        return
    logging.info("Slicing BAM %s:%d-%d → %s", chrom, start, end, out_bam)
    with open(out_bam, "wb") as fh:
        run(["samtools", "view", "-b", str(full_bam), f"{chrom}:{start}-{end}"], stdout=fh)
    run(["samtools", "index", str(out_bam)])


def awk_slice(in_path: Path, chrom: str, start: int, end: int,
              out_path: Path, col_start: int = 2, col_end: int = 3,
              dry: bool = False):
    """Slice text-based BED/TSV by region using awk."""
    if dry:
        logging.info("[DRY] awk slice %s on %s:%d-%d → %s",
                     in_path.name, chrom, start, end, out_path)
        return
    if out_path.exists():
        logging.info("Slice exists (%s), skipping", out_path)
        return
    logging.info("Slicing %s on %s:%d-%d → %s",
                 in_path.name, chrom, start, end, out_path)
    cmd = (
        f"awk '$1==\"{chrom}\" && ${col_start}>={start} && ${col_end}<={end}' "
        f"{shlex.quote(str(in_path))} > {shlex.quote(str(out_path))}"
    )
    run(["bash", "-c", cmd])


def slice_gtf_fully_contained(gtf_path: Path, chrom: str, start: int,
                              end: int, out_path: Path,
                              dry: bool = False):
    """Slice GTF features fully contained in the region."""
    if dry:
        logging.info("[DRY] slice GTF %s on %s:%d-%d → %s",
                     gtf_path.name, chrom, start, end, out_path)
        return
    if out_path.exists():
        logging.info("GTF slice exists (%s), skipping", out_path)
        return
    logging.info("Slicing GTF %s on %s:%d-%d → %s",
                 gtf_path.name, chrom, start, end, out_path)
    import re
    target = chrom if chrom.startswith("chr") else f"chr{chrom}"
    feats, lines, genes = {}, {}, []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] != target:
                continue
            typ, st, en = cols[2], int(cols[3]), int(cols[4])
            if typ == "gene":
                genes.append((st, en, line))
            m = re.search(r'transcript_id\s+"([^"]+)"', cols[8])
            if not m:
                continue
            tid = m.group(1)
            feats.setdefault(tid, []).append((st, en))
            lines.setdefault(tid, []).append(line)
    overlaps = [
        iv
        for ivs in feats.values()
        for iv in ivs
        if not (iv[1] < start or iv[0] > end)
    ]
    ns, ne = ((min(iv[0] for iv in overlaps), max(iv[1] for iv in overlaps))
              if overlaps else (start, end))
    with open(out_path, "w") as out:
        for g_st, g_en, g_line in sorted(genes):
            if g_st >= ns and g_en <= ne:
                out.write(g_line)
        for tid, ivs in feats.items():
            if all(iv[0] >= ns and iv[1] <= ne for iv in ivs):
                for line in sorted(lines[tid], key=lambda l: int(l.split("\t")[3])):
                    out.write(line)


def mk_region_fasta(bed_path: Path, fa_full: Path, fa_out: Path, dry: bool = False):
    """Build a FASTA slice of reads overlapping the BED regions."""
    if dry:
        logging.info("[DRY] mk_region_fasta from %s → %s", fa_full.name, fa_out.name)
        return
    if fa_out.exists():
        logging.info("FASTA exists (%s), skipping", fa_out)
        return
    logging.info("Building FASTA slice from %s → %s", fa_full.name, fa_out.name)
    ids = set()
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                ids.add(parts[3])
    write = False
    with open(fa_full) as inf, open(fa_out, "w") as outf:
        for line in inf:
            if line.startswith(">"):
                hdr = line[1:].split()[0]
                write = hdr in ids
            if write:
                outf.write(line)


def slice_region(full_bam: Path, query_bed: Path, gtf: Path, junctions,
                 reads, exp5_bed, exp3_bed,
                 ref5_bed, ref3_bed,
                 chrom: str, start: int, end: int,
                 outdir: Path, dry: bool = False):
    """Slice all inputs for one region into `outdir`."""
    outdir.mkdir(parents=True, exist_ok=True)
    logging.info("→ slice_region into %s", outdir)

    # 1) BAM
    slice_bam(full_bam, chrom, start, end, outdir / f"{chrom}_{start}-{end}.bam", dry)

    # 2) query BED
    awk_slice(query_bed, chrom, start, end, outdir / f"{chrom}_{start}-{end}.bed", dry=dry)

    # 3) experimental beds
    if exp5_bed:
        awk_slice(exp5_bed, chrom, start, end, outdir / f"{chrom}_{start}-{end}.exp5.bed", dry=dry)
    if exp3_bed:
        awk_slice(exp3_bed, chrom, start, end, outdir / f"{chrom}_{start}-{end}.exp3.bed", dry=dry)

    # 4) reference beds
    if ref5_bed:
        awk_slice(ref5_bed, chrom, start, end, outdir / f"{chrom}_{start}-{end}.ref5.bed", dry=dry)
    if ref3_bed:
        awk_slice(ref3_bed, chrom, start, end, outdir / f"{chrom}_{start}-{end}.ref3.bed", dry=dry)

    # 5) GTF
    slice_gtf_fully_contained(gtf, chrom, start, end, outdir / f"{chrom}_{start}-{end}.gtf", dry=dry)

    # 6) junctions
    if junctions:
        awk_slice(Path(junctions), chrom, start, end, outdir / f"{chrom}_{start}-{end}.tab", dry=dry)

    # 7) FASTA reads
    fasta_list = reads if isinstance(reads, (list, tuple)) else [reads]
    for fa in fasta_list:
        mk_region_fasta(outdir / f"{chrom}_{start}-{end}.bed",
                        Path(fa),
                        outdir / f"{chrom}_{start}-{end}.fasta",
                        dry=dry)


def slice_all_regions(cfg, dry_run: bool = False):
    """Slice all configured regions for each corrected run under `outputs/<sample_name>/…`."""
    # Only run if not whole-sample and regions defined
    if getattr(cfg, 'whole_sample', False):
        logging.info("Whole-sample mode: skipping region slicing")
        return
    regions = getattr(cfg.run, 'regions', [])
    if not regions:
        logging.info("No regions configured: skipping region slicing")
        return

    # Base output dir from config
    base_out    = Path("outputs") / cfg.sample_name
    correct_base = base_out / "correct"

    # Resolve inputs
    input_root  = Path(cfg.input_root)
    data        = cfg.data_dir
    gtf         = input_root / data.annotation_gtf
    junctions   = (input_root / data.junctions) if getattr(data, 'junctions', None) else None
    reads_raw   = data.reads_fasta
    reads       = ([input_root / reads_raw]
                   if isinstance(reads_raw, str)
                   else [input_root / rf for rf in reads_raw])
    exp5_bed    = (input_root / data.experimental_tss_bed) if getattr(data, 'experimental_tss_bed', None) else None
    exp3_bed    = (input_root / data.experimental_tts_bed) if getattr(data, 'experimental_tts_bed', None) else None
    ref5_bed    = (input_root / data.ref_tss_bed) if getattr(data, 'ref_tss_bed', None) else None
    ref3_bed    = (input_root / data.ref_tts_bed) if getattr(data, 'ref_tts_bed', None) else None

    # Look for each corrected run
    pattern     = str(correct_base / "*" / "*")
    runs = sorted(glob.glob(pattern))
    if not runs:
        msg = f"No corrected runs found under {correct_base}; cannot slice regions."
        logging.error(msg)
        raise RuntimeError(msg)
    for sample_run in sorted(glob.glob(pattern)):
        sample   = Path(sample_run).parent.name
        corr_run = Path(sample_run).name

        # find the corrected BAM/BED
        bam_files = glob.glob(str(Path(sample_run) / f"{sample}*.bam"))
        bed_files = glob.glob(str(Path(sample_run) / f"{sample}*.bed"))
        if not bam_files or not bed_files:
            logging.warning("No BAM/BED in %s, skipping", sample_run)
            continue

        full_bam  = Path(bam_files[0])
        query_bed = Path(bed_files[0])

        # now slice each region
        for region in regions:
            chr_ = region.chr
            for span in region.ranges:
                start, end = map(int, span.split(":", 1))
                region_id  = f"{chr_}_{start}_{end}"
                outdir_r   = base_out / "regions" / region_id / sample / corr_run
                slice_region(full_bam, query_bed, gtf, junctions,
                             reads, exp5_bed, exp3_bed,
                             ref5_bed, ref3_bed,
                             chr_, start, end, outdir_r, dry_run)
