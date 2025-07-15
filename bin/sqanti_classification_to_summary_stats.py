#!/usr/bin/env python3
import sys
import os
import matplotlib.pyplot as plt
from math import isnan

"""
Usage:
    python sqanti3_classification_summary.py <manifest.txt> <summary_prefix>

<manifest.txt> lists one path to a *_classification.txt per line.
"""

# --- parse args -----------------------------------------------------------
if len(sys.argv) != 3:
    sys.exit("Usage: python sqanti3_classification_summary.py <manifest.txt> <summary_prefix>")
manifestfile = sys.argv[1]
summaryprefix = sys.argv[2]

# --- 1) Read file list, skip Mac resource forks & non-existent files -------
files = []
with open(manifestfile, 'r') as mf:
    for l in mf:
        path = l.strip()
        if not path:
            continue
        bn = os.path.basename(path)
        if bn.startswith('._'):
            continue
        if not os.path.isfile(path):
            print(f"[WARN] file not found, skipping: {path}")
            continue
        files.append(path)

if not files:
    sys.exit("[ERROR] No valid classification files found in manifest.")

# --- 2) Canonical SQANTI classes -----------------------------------------
isoclasses = [
    'FSM', 'ISM', 'NIC', 'NNC',
    'Genic_Genomic', 'Genic_Intron',
    'Antisense', 'Fusion', 'Intergenic'
]

# --- 3) Raw→canonical mapping ---------------------------------------------
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
normalized_raw_to_sqanti = {
    k.strip().lower().replace('-', '_'): v
    for k, v in _raw_keys.items()
}

# --- 4) Helpers ------------------------------------------------------------
def normalize_label(lbl):
    return lbl.strip().lower().replace('-', '_')

def safe_float(v):
    try:
        return float(v)
    except:
        return float('nan')

# --- 5) Main loop: compute summary per file -------------------------------
alllines = []

for path in files:
    tot = cages = polya = ref5 = ref3 = srtm = sntm = 0
    counts = {c: 0 for c in isoclasses}
    warned = set()

    with open(path, encoding='utf-8', errors='replace') as f:
        for row in f:
            cols = row.rstrip('\n').split('\t')
            if cols[0] == 'isoform' or len(cols) < 44:
                continue

            tot += 1
            raw = cols[5]
            norm = normalize_label(raw)
            iso = normalized_raw_to_sqanti.get(norm)
            if iso in counts:
                counts[iso] += 1
            else:
                if norm not in warned:
                    print(f"[WARN] Unknown class '{raw}' in {path}")
                    warned.add(norm)

            # CAGE/polyA support
            dC = safe_float(cols[39])
            iC = cols[40].lower() == 'true'
            dP = safe_float(cols[41])
            iP = cols[42].lower() == 'true'
            hasC = iC or (not isnan(dC) and abs(dC) <= 50)
            hasP = iP or (not isnan(dP) and abs(dP) <= 50)

            # Annotated TSS/TTS support
            dTSS = safe_float(cols[10])
            dTTS = safe_float(cols[11])
            gTSS = safe_float(cols[12])
            gTTS = safe_float(cols[13])
            a5 = ((not isnan(dTSS) and abs(dTSS) <= 50) or (not isnan(gTSS) and abs(gTSS) <= 50))
            a3 = ((not isnan(dTTS) and abs(dTTS) <= 50) or (not isnan(gTTS) and abs(gTTS) <= 50))

            if a5:
                ref5 += 1
            if a3:
                ref3 += 1

            if iso in ('FSM', 'ISM'):
                if (hasC or a5) and (hasP or a3):
                    srtm += 1
            else:
                if (hasC or a5) and (hasP or a3):
                    sntm += 1

            cages += hasC
            polya += hasP

    summary = [
        path.rsplit('_classification.txt', 1)[0],
        round(srtm / tot, 3) if tot else 0,
        round(sntm / tot, 3) if tot else 0,
        round(cages / tot, 3) if tot else 0,
        round(polya / tot, 3) if tot else 0,
        round(ref5 / tot, 3) if tot else 0,
        round(ref3 / tot, 3) if tot else 0
    ] + [counts[c] for c in isoclasses] + [tot]
    alllines.append(summary)

# --- 6) Write TSV ---------------------------------------------------------
header = [
    'filename','fracSRTM','fracSNTM',
    "5'endsupport","3'endsupport",
    "5'refsupport","3'refsupport"
] + isoclasses + ['total_isoforms']
with open(f"{summaryprefix}.tsv", 'w') as out:
    out.write('\t'.join(header) + '\n')
    for line in alllines:
        out.write('\t'.join(map(str, line)) + '\n')

# --- 7) Plot results ------------------------------------------------------
num_files = len(alllines)
fig_width = max(20, num_files * 0.2)  # 0.2 inches per file, minimum 20 inches
fig, axs = plt.subplots(4, 1, figsize=(fig_width, 13))
x = list(range(num_files))

# Panel 1: SRTM/SNTM
s1 = [l[1] for l in alllines]
s2 = [l[2] for l in alllines]
axs[0].bar(x, s1, label='SRTM')
axs[0].bar(x, s2, bottom=s1, label='SNTM')
axs[0].set_ylim(0, 1)
axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[0].set_xticks(x, [''] * num_files)
axs[0].set_ylabel("SRTM/SNTM")

# Panel 2: 5'/3' end support (CAGE/polyA)
c5 = [l[3] for l in alllines]
c3 = [l[4] for l in alllines]
axs[1].bar([i - 0.2 for i in x], c5, width=0.4, label="5'endsupport")
axs[1].bar([i + 0.2 for i in x], c3, width=0.4, label="3'endsupport")
axs[1].set_ylim(0, 1)
axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[1].set_xticks(x, [''] * num_files)
axs[1].set_ylabel("CAGE/polyA")

# Panel 3: 5'/3' reference support
r5 = [l[5] for l in alllines]
r3 = [l[6] for l in alllines]
axs[2].bar([i - 0.2 for i in x], r5, width=0.4, label="5'refsupport")
axs[2].bar([i + 0.2 for i in x], r3, width=0.4, label="3'refsupport")
axs[2].set_ylim(0, 1)
axs[2].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[2].set_xticks(x, [''] * num_files)
axs[2].set_ylabel("Ref TSS/TTS")

# Panel 4: isoform class counts
bottom = [0] * num_files
for idx, cls in enumerate(isoclasses):
    vals = [l[7 + idx] for l in alllines]
    axs[3].bar(x, vals, bottom=bottom, label=cls)
    bottom = [b + v for b, v in zip(bottom, vals)]
axs[3].set_xticks(x)
axs[3].set_xticklabels(
    [l[0].split('/')[-1] for l in alllines],
    rotation=90, ha='center', fontsize=6
)
axs[3].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[3].set_ylabel("Isoform count")

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig(f"{summaryprefix}.png", dpi=600)



'''

python /private/groups/brookslab/hdheath/projects/rotation/test_suite/sqanti_classification_to_summary_stats_v.4.py \
     aim-1-for-sqanti-summary \
     sqanti_aim-1

python /private/groups/brookslab/hdheath/projects/rotation/test_suite/sqanti_classification_to_summary_stats_v.4.py \
     yeast-for-sqanti \
     sqanti_aim-2_yeast

python /private/groups/brookslab/hdheath/projects/rotation/test_suite/sqanti_classification_to_summary_stats_v.4.py \
     human-for-sqanti_summary-stats \
     sqanti_aim-2_human


'''