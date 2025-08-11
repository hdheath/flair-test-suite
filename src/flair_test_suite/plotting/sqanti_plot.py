#!/usr/bin/env python3
"""
sqanti_plot.py

Find every sqanti_results.tsv under results/, plot a 3‐panel figure,
and write sqanti.png under plots/region/run_name/, skipping if present.
"""

import os, sys, glob, argparse, pandas as pd, matplotlib.pyplot as plt

# ─── same hack ───────────────────────────────────────────────────────────────
THIS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(THIS, "..", ".."))
SRC  = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

def plot_summary(summary_tsv, plot_dir):
    out_png = os.path.join(plot_dir, 'sqanti.png')
    if os.path.exists(out_png):
        print(f"[SKIP] plot already exists: {out_png}")
        return

    os.makedirs(plot_dir, exist_ok=True)
    df = pd.read_csv(summary_tsv, sep='\t')
    samples = df['sample'].tolist()
    x       = range(len(df))

    # panel 1
    fig, axs = plt.subplots(3,1, figsize=(max(8,len(df)*0.4), 10), constrained_layout=True)
    axs[0].bar(x, df['fracSRTM'], label='SRTM')
    axs[0].bar(x, df['fracSNTM'], bottom=df['fracSRTM'], label='SNTM')
    axs[0].set_ylabel("SRTM/SNTM"); axs[0].legend(loc='center left', bbox_to_anchor=(1,0.5)); axs[0].set_xticks([])

    # panel 2
    axs[1].bar([i-0.2 for i in x], df["5'endsupport"], width=0.4, label="5'")
    axs[1].bar([i+0.2 for i in x], df["3'endsupport"], width=0.4, label="3'")
    axs[1].set_ylabel("end support"); axs[1].legend(loc='center left', bbox_to_anchor=(1,0.5)); axs[1].set_xticks([])

    # panel 3: isoform counts
    class_cols = [c for c in df.columns if c in ISOCLASSES]
    bottom = [0]*len(df)
    for cls in class_cols:
        axs[2].bar(x, df[cls], bottom=bottom, label=cls)
        bottom = [b+v for b,v in zip(bottom, df[cls])]
    axs[2].set_xticks(x)
    axs[2].set_xticklabels(samples, rotation=90, fontsize=6)
    axs[2].set_ylabel("isoform count")
    axs[2].legend(loc='center left', bbox_to_anchor=(1,0.5))

    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"[✓] Wrote plot {out_png}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-o','--outdir', required=True)
    args = p.parse_args()

    pattern = os.path.join(args.outdir,'results','*','*','sqanti_results.tsv')
    for tsv in sorted(glob.glob(pattern)):
        region = os.path.basename(os.path.dirname(os.path.dirname(tsv)))
        runn   = os.path.basename(os.path.dirname(tsv))
        plot_dir = os.path.join(args.outdir,'plots',region,runn)
        print(f"Plotting {tsv} → {plot_dir}")
        plot_summary(tsv, plot_dir)

if __name__ == '__main__':
    main()
