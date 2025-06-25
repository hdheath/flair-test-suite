# FLAIR Test Suite

A testing pipeline to:

- **FLAIR Align** with user inputted options
- **FLAIR Correct** with user inputted options
- **Slice** input data into user inputted region-specific files  
- **FLAIR Collapse** each region with user inputted options  
- **FLAIR Collapse QC** with SQANTI3 and **plot** results  
- **FLAIR Collapse QC** with TED and **plot** results  
- **FLAIR Quantify** each collapse run with user inputted options  
- **FLAIR Quantify QC** 

---

## Pipeline Workflow

### **FLAIR Align**
**Note - this is only ran for version <3**

**Inputs:**
- `reference.fasta`
- `reads.fasta`

**Outputs:**
- `aligned.bed`
- `aligned.bam`
- `reads.fasta`

---

### **FLAIR Correct**
**Note - this is only ran for version <3**

**Inputs:**
- `aligned.bed`
- `aligned.bam`
- `reference.fasta`
- *optional*: `reference.gtf`
- *optional*: `reference_splice_junctions.txt`

**Outputs:**
- `corrected_aligned.bed`
- `corrected_aligned.bam`

---

### **Slice**

**Inputs:**
- Region(s) of interest: `chr:start-end`
- `reads.fasta`
- `corrected_aligned.bed`
- `corrected_aligned.bam`
- *optional*: `reference.gtf`
- *optional*: TSS/TTS BED files (reference or experimental)

**Outputs:**
- Sliced files of each region from the input files  
  **Note:** Only features *fully contained* within the given region of interest will be included in the output file

---

### **FLAIR Collapse**

**Inputs:**
- `reference.fasta`
- Sliced `reference.gtf`
- Sliced, corrected, aligned `bed` file
- Sliced, corrected, aligned `reads.fasta` file

**Outputs:**
- (TBD)

---

### **Collapse QC – SQANTI**

**Inputs:**
- (TBD)

**Outputs:**
- (TBD)

---

### **Collapse QC – TED**

**Inputs:**
- (TBD)

**Outputs:**
- (TBD)

---

### **FLAIR Quantify**

**Inputs:**
- (TBD)

**Outputs:**
- (TBD)

---

### **Quantify QC**

**Inputs:**
- (TBD)

**Outputs:**
- (TBD)



## Pipeline Workflow Schematic 


## Example Directory Layout
```plaintext
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
│               ├── sqanti_error.log
│               └── quantify_qc.log
└── tests/
    └── data/
        ├── sample.gtf
        └── sample_reads.fastq
        └── reference.fastq
        └── [optional] reference TSS/TTS bed file(s)
        └── [optional] experiment TSS/TTS bed file(s)
        └── [optional] reference splice junction bed file
```


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

