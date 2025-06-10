#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
import sys
import matplotlib.pyplot as plt
from math import isnan

def run_sqanti_qc(gtf, ref_gtf, ref_fasta, prime5, prime3, sj, prefix, out_subdir, cpus, dry_run):
    """Invoke SQANTI3 QC for one sample."""
    cmd = [
        "sqanti3_qc",
        gtf,
        ref_gtf,
        ref_fasta,
        "--CAGE_peak", prime5,
        "--polyA_peak", prime3,
        "-c", sj,
        "-o", prefix,
        "-d", out_subdir,
        "--skipORF",
        "--cpus", str(cpus)
    ]
    print(f"{'[DRY RUN]' if dry_run else '[RUNNING]'} {' '.join(cmd)}")
    if not dry_run:
        subprocess.run(cmd, check=True)

def parse_args():
    p = argparse.ArgumentParser(
        description="Run SQANTI3 QC on a manifest then summarize all classification.txt files."
    )
    p.add_argument("-m","--manifest", required=True,
                   help="Manifest TSV: columns gtf, ref_gtf, ref_fasta, prime5, prime3, sj")
    p.add_argument("-d","--outdir", required=True,
                   help="Base dir for per-sample SQANTI runs")
    p.add_argument("-c","--cpus", type=int, default=8,
                   help="CPUs to pass to SQANTI3")
    p.add_argument("-s","--summary_prefix", required=True,
                   help="Prefix (no extension) for summary .tsv and .png")
    p.add_argument("--dry_run", action="store_true",
                   help="Only print commands, don’t actually run SQANTI3")
    return p.parse_args()

# canonical classes and mapping
ISOCLASSES = [
    'FSM', 'ISM', 'NIC', 'NNC',
    'Genic_Genomic', 'Genic_Intron',
    'Antisense', 'Fusion', 'Intergenic'
]
_raw_keys = {
    'full-splice_match':       'FSM',
    'incomplete-splice_match': 'ISM',
    'novel_in_catalog':        'NIC',
    'novel-not-in-catalog':    'NNC',
    'novel_not_in_catalog':    'NNC',
    'genic':                   'Genic_Genomic',
    'genic_intron':            'Genic_Intron',
    'genic_genomic':           'Genic_Genomic',
    'antisense':               'Antisense',
    'fusion':                  'Fusion',
    'intergenic':              'Intergenic'
}
NORMALIZED_MAP = {
    k.strip().lower().replace('-', '_'): v
    for k, v in _raw_keys.items()
}

def normalize_label(lbl):
    return lbl.strip().lower().replace('-', '_')

def safe_float(x):
    try:
        return float(x)
    except:
        return float('nan')

def summarize_classifications(class_files, summary_prefix):
    files = []
    for path in class_files:
        bn = os.path.basename(path)
        if bn.startswith('._'): continue
        if not os.path.isfile(path):
            print(f"[WARN] missing, skipping: {path}", file=sys.stderr)
            continue
        files.append(path)
    if not files:
        sys.exit("[ERROR] No classification files found to summarize.")

    all_summaries = []
    for path in files:
        tot = cages = polya = srtm = sntm = 0
        counts = {c:0 for c in ISOCLASSES}
        warned = set()

        with open(path, encoding='utf-8', errors='replace') as fh:
            for row in fh:
                cols = row.rstrip('\n').split('\t')
                if cols[0]=='isoform' or len(cols)<44:
                    continue
                tot += 1

                # class
                raw = cols[5]
                iso = NORMALIZED_MAP.get(normalize_label(raw))
                if iso in counts:
                    counts[iso] += 1
                else:
                    if normalize_label(raw) not in warned:
                        print(f"[WARN] Unknown class '{raw}' in {path}", file=sys.stderr)
                        warned.add(normalize_label(raw))

                # 5'/3' support
                dC = safe_float(cols[39]); iC = cols[40].lower()=='true'
                dP = safe_float(cols[41]); iP = cols[42].lower()=='true'
                hasC = iC or (not isnan(dC) and abs(dC)<=50)
                hasP = iP or (not isnan(dP) and abs(dP)<=50)

                # annotated TSS/TTS support
                dTSS = safe_float(cols[10]); dTTS = safe_float(cols[11])
                gTSS = safe_float(cols[12]); gTTS = safe_float(cols[13])
                a5 = (not isnan(dTSS) and abs(dTSS)<=50) or (not isnan(dTTS) and abs(dTTS)<=50)
                a3 = (not isnan(gTSS) and abs(gTSS)<=50) or (not isnan(gTTS) and abs(gTTS)<=50)

                # SRTM vs SNTM
                if iso in ('FSM','ISM'):
                    if (hasC or a5) and (hasP or a3):
                        srtm += 1
                else:
                    if (hasC or a5) and (hasP or a3):
                        sntm += 1

                cages += hasC
                polya += hasP

        base = path.rsplit('_classification.txt',1)[0]
        summary = [
            base,
            round(srtm/tot,3) if tot else 0,
            round(sntm/tot,3) if tot else 0,
            round(cages/tot,3) if tot else 0,
            round(polya/tot,3) if tot else 0
        ] + [counts[c] for c in ISOCLASSES] + [tot]
        all_summaries.append(summary)

    # write TSV
    header = ['filename','fracSRTM','fracSNTM',"5'endsupport","3'endsupport"] + ISOCLASSES + ['total_isoforms']
    with open(f"{summary_prefix}.tsv",'w') as out:
        out.write('\t'.join(header)+'\n')
        for row in all_summaries:
            out.write('\t'.join(map(str,row))+'\n')

    # plot
    fig, axs = plt.subplots(3,1, figsize=(6,10))
    x = range(len(all_summaries))

    # panel 1
    s1 = [r[1] for r in all_summaries]
    s2 = [r[2] for r in all_summaries]
    axs[0].bar(x, s1, label='SRTM')
    axs[0].bar(x, s2, bottom=s1, label='SNTM')
    axs[0].set_ylim(0,1)
    axs[0].legend(loc='center left', bbox_to_anchor=(1,0.5))
    axs[0].set_xticks(x, ['']*len(x))

    # panel 2
    c5 = [r[3] for r in all_summaries]
    c3 = [r[4] for r in all_summaries]
    axs[1].bar([i-0.2 for i in x], c5, width=0.3, label="5' support")
    axs[1].bar([i+0.2 for i in x], c3, width=0.3, label="3' support")
    axs[1].set_ylim(0,1)
    axs[1].legend(loc='center left', bbox_to_anchor=(1,0.5))
    axs[1].set_xticks(x, ['']*len(x))

    # panel 3
    bottom = [0]*len(x)
    for idx, cls in enumerate(ISOCLASSES):
        vals = [r[5+idx] for r in all_summaries]
        axs[2].bar(x, vals, bottom=bottom, label=cls)
        bottom = [b+v for b,v in zip(bottom, vals)]
    axs[2].set_xticks(x, [os.path.basename(r[0]) for r in all_summaries], rotation=90, ha='center')
    axs[2].legend(loc='center left', bbox_to_anchor=(1,0.5))

    plt.tight_layout()
    plt.savefig(f"{summary_prefix}.png", dpi=600)
    print(f"[✓] Summary written to {summary_prefix}.tsv and {summary_prefix}.png")

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    classification_files = []
    # 1) Run SQANTI3 QC
    with open(args.manifest, newline='') as mf:
        reader = csv.DictReader(mf, delimiter='\t')
        for row in reader:
            gtf       = row["gtf"]
            ref_gtf   = row["ref_gtf"]
            ref_fasta = row["ref_fasta"]
            prime5    = row["prime5"]
            prime3    = row["prime3"]
            sj        = row["sj"]

            prefix = os.path.basename(gtf).replace(".isoforms.gtf","")
            out_subdir = os.path.join(args.outdir, prefix)
            os.makedirs(out_subdir, exist_ok=True)

            run_sqanti_qc(
                gtf, ref_gtf, ref_fasta, prime5, prime3, sj,
                prefix, out_subdir, args.cpus, args.dry_run
            )

            # record where the classification file will be
            classification_files.append(
                os.path.join(out_subdir, f"{prefix}_classification.txt")
            )

    # 2) Summarize all classification.txt files
    summarize_classifications(classification_files, args.summary_prefix)

if __name__ == "__main__":
    main()




'''

python automated-sqanti-runs_v.4.py \
  -m yeast-for-sqanti.tsv \
  -d /private/groups/brookslab/hdheath/projects/rotation/test_suite/sqanti/aim-2/yeast \
  -c 12 \
  -s yeast-aim-2_sqanti_summary \
  --dry_run


python automated-sqanti-runs_v.4.py \
  -m human-for-sqanti.tsv \
  -d /private/groups/brookslab/hdheath/projects/rotation/test_suite/sqanti/aim-2/human \
  -c 12 \
  -s human-aim-2_sqanti_summary \
  --dry_run



'''