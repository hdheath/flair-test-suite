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
│   └── manifest.py            # pipeline parameters as a Python dict
├── src/
│   └── flair_automate/        # core library
│       ├── __init__.py
│       ├── flair_align_automation.py
│       ├── flair_correct_automation.py
│       ├── slicer.py
│       ├── flair_collapse_automation.py
│       ├── parse_region_metrics_from_gtf.py
│       ├── sqanti_runner.py
│       └── sqanti_plot.py
├── scripts/
│   └── flair_test_suite.py    # top-level “run everything” CLI
├── outputs/                   # all generated data
│   ├── align/                 # FLAIR-align outputs
│   │   └── <sample>/
│   │       └── align1/        # or “default” if no flags, or a name derived from flags
│   │           ├── <sample>.bam
│   │           └── <sample>.bed
│   ├── correct/               # FLAIR-correct outputs
│   │   └── <sample>/
│   │       └── <align_run>/   
│   │           └── <corr_run>/     # or name derived from flags
│   │               ├── <sample>.corrected.bam
│   │               └── <sample>.corrected.bed
│   ├── regions/               # per-region slices for every sample/run
│   │   └── chr20_3218000_3250000/
│   │       └── <sample>/
│   │           └── <corr_run>/  
│   │               ├── chr20_3218000-3250000.bam
│   │               ├── chr20_3218000-3250000.bed
│   │               ├── chr20_3218000-3250000.gtf
│   │               ├── chr20_3218000-3250000.fasta
│   │               ├── chr20_3218000-3250000.exp5.bed
│   │               ├── chr20_3218000-3250000.exp3.bed
│   │               ├── chr20_3218000-3250000.ref5.bed
│   │               └── chr20_3218000-3250000.ref3.bed
│   ├── results/               # FLAIR-collapse isoforms
│   │   └── chr20_3218000_3250000/
│   │       └── <sample>/
│   │           └── <corr_run>/
│   │               ├── <collapse_run>/                       # or named by collapse flags
│   │               │   ├── *.isoforms.gtf
│   │               │   └── *.isoforms.bed
│   │               └── <collapse_run>/ # example flag-derived name
│   │                   ├── *.isoforms.gtf
│   │                   └── *.isoforms.bed
│   ├── logs/                  # timing & error logs
│   │   ├── chr20_3218000_3250000/
│   │   │   └── <sample>/
│   │   │       └── <corr_run>/
                    correct.time.log
│   │   │           └── <collapse_run>/              # e.g. “default”
│   │   │               └── collapse.time.log
│   │   └── sqanti/            # SQANTI errors
│   │       └── chr20_3218000_3250000/
│   │           └── <sample>/
│   │               └── <corr_run>/
│   │                 └── <collapse_run>/
│   │                       └── sqanti_error.log
│   └── plots/                 # SQANTI QC figures
│       └── chr20_3218000_3250000/
│           └── <sample>/
│               └── <corr_run>/
│                    └── <collapse_run>/
│                             └── region_metrics.png
│                             └── sqanti.png
│                             └── sqanti_summary.tsv
│                             └── ted.png
│                             └── ted_summary.tsv
└── tests/
    └── data/                  # data for unit tests
        ├── sample.gtf
        └── sample_reads.fastq




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

