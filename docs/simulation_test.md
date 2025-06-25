# FLAIR Test Suite

A testing pipeline to:

- **FLAIR Align** with user inputted options
- **FLAIR Correct** with user inputted options
- **Slice** input data into user inputted region-specific files  
- **FLAIR Collapse** each region with user inputted options  
- **QC** with SQANTI3 and **plot** results  
- **QC** with TED and **plot** results  

---

## Steps

**FLAIR Align**

**FLAIR Correct**

**Slice**

**FLAIR Collapse**

**QC**


## Example Directory Layout

flair-test-suite/
├── config/
│   └── manifest.py
├── src/
│   └── flair_automate/
│       ├── __init__.py
│       ├── flair_align_automation.py
│       ├── flair_correct_automation.py
│       ├── slicer.py
│       ├── flair_collapse_automation.py
│       ├── parse_region_metrics_from_gtf.py
│       ├── sqanti_runner.py
│       └── sqanti_plot.py
├── scripts/
│   └── flair_test_suite.py
├── outputs/
│   └── <sample>/                                # organism, fasta, version
│       └── run_<align_id>[_<corr_id>]/          # Alignment settings, Correct settings
│           ├── align/                                 # FLAIR align output 
│           │   ├── <sample>.bam
│           │   └── <sample>.bed
│           │
│           ├── correct/                               # FLAIR correct output
│           │   ├── <sample>.corrected.bam
│           │   └── <sample>.corrected.bed
│           │
│           ├── regions/                               # sliced files used for FLAIR collapse and FLAIR quantify runs
│           │   └── chr20_3218000_3250000/
│           │       ├── raw/                               
│           │       │   ├── chr20…-3250000.bam
│           │       │   ├── chr20…-3250000.bed
│           │       │   ├── chr20…-3250000.gtf
│           │       │   ├── chr20…-3250000.fasta
│           │       │   ├── chr20…-3250000.exp5.bed
│           │       │   └── …  
│           │       │
│           │       ├── collapse/                           # FLAIR collapse per-region, using different options 
│           │       │   └── collapse_<flags>/
│           │       │       ├── *.isoforms.gtf
│           │       │       └── *.isoforms.bed
│           │       │
│           │       ├── collapse_qc/                        # FLAIR collapse QC for this region across all runs 
│           │       │   ├── sqanti_summary.tsv
│           │       │   ├── sqanti.png
│           │       │   ├── ted_summary.tsv
│           │       │   ├── ted.png
│           │       │   └── region_metrics.png
│           │       │
│           │       └── quantify_<flags>/                   # FLAIR quantify per-region, using different options 
│           │           ├── isoform_tpm.tsv
│           │           └── gene_counts.tsv
│           │       └── quantify_qc/                        # FLAIR quantify QC for this region across all runs 
│           │           ├── flair_quantify_metrics.png
│           │           
│           │
│           └── logs/
│               ├── align.time.log
│               ├── correct.time.log
│               ├── collapse.time.log
│               ├── quantify.time.log
│               ├── ted.time.log
│               └── sqanti_error.log
└── tests/
    └── data/
        ├── sample.gtf
        └── sample_reads.fastq
        └── reference.fastq
        └── [optional] reference TSS/TTS bed file(s)
        └── [optional] experiment TSS/TTS bed file(s)
        └── [optional] reference splice junction bed file






## Pipeline Workflow Schematic 
**convert this to bullet lists** 
                ┌─────────────────────────────────┐
                │      scripts/driver.py         │
                └──────────────┬──────────────────┘
                               ▼
╔═══════════════════════╗    1) slice_all_regions(cfg)   ╔══════════════════════╗
║  config/manifest.yaml ║ ──────────────────────────────▶║ outputs/regions/      ║
╚═══════════════════════╝                                 ║ └─ chr20_3218000_3250000/  ║
                                                          ║     ├─ chr20_3218000-3250000.bam   ║
                                                          ║     ├─ chr20_3218000-3250000.bed   ║
                                                          ║     ├─ … (.gtf, .fasta, .5prime.bed)║
                                                          ╚══════════════════════════════════╝
                                                               │
                                                               ▼
                                                    2) parse_all_regions(…)
                                                               │
                                                          ╔═════════════════════╗
                                                          ║ outputs/regions/…   ║
                                                          ║ └─ chr…/            ║
                                                          ║     ├─ gene_summary.csv      ║
                                                          ║     └─ transcript_summary.csv║
                                                          ╚═════════════════════╝
                                                               │
                                                               ▼
                                                    3) collapse_all_regions(cfg)
                                                               │
                                                          ╔════════════════════════╗
                                                          ║ outputs/results/       ║
                                                          ║ └─ chr20_…_/           ║
                                                          ║     ├─ default/        ║
                                                          ║     │   └─ *.isoforms.gtf   ║
                                                          ║     └─ <flag_combo>/   ║
                                                          ╚════════════════════════╝
                                                               │
                                                               ▼
                                       (Optional) 4) SQANTI QC & summary
                                                               │
                                                          ╔══════════════════════╗
                                                          ║ outputs/results/…    ║
                                                          ║ └─ chr…_/            ║
                                                          ║     └─ <run>/        ║
                                                          ║         ├─ *_classification.txt  ║
                                                          ║         └─ sqanti_results.tsv    ║
                                                          ║ outputs/logs/sqanti/              ║
                                                          ║ └─ <run>/                       ║
                                                          ║     ├─ sqanti.log               ║
                                                          ║     └─ errors.log               ║
                                                          ╚══════════════════════╝
                                                               │
                                                               ▼
                                       (Optional) 5) SQANTI Plotting
                                                               │
                                                          ╔══════════════════════╗
                                                          ║ outputs/plots/       ║
                                                          ║ └─ chr…_/            ║
                                                          ║     └─ <run>/        ║
                                                          ║         └─ sqanti.png       ║
                                                          ╚══════════════════════╝

## Extending & Testing

* add analysis to pipeline
   * need to downlaod sqanti to flair env : run sqanti w/ per region plot 
      * save to results/region/run_name/sqanti_results.tsv and plots/region/run_name/sqanti.png
   * TED w/ per region plot - compare with gtf metrics - parallel coordinate map on top of box plots
      * save to results/region/run_name/ted_results.tsv and plots/region/run_name/ted.png
* make fasta reads generator faster (wont need this with new version)
* make testing for flair dev (new flair automation py)
* Add new regions to `config/manifest.yaml`.
* Update `tests/data/` fixtures for unit/integration tests.
* Run `pytest` to verify functions and end-to-end behavior.
* create optional plotting of regions

