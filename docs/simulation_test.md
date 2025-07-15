# FLAIR Test Suite Workflow

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
**Note - this is only ran for version < 3 - newer FLAIR versions allow users to choose their alignment method**

**Inputs:**
- `reference.fasta`
- `reads.fasta`

**Outputs:**
- `aligned.bed`
- `aligned.bam`
- `reads.fasta`

---

### **FLAIR Correct**
**Note - this is only ran for version < 3, correct is within newer FLAIR versions' collapse command**

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

**Inputs for version < 3**
- `reference.fasta`
- Sliced `reference.gtf`
- Sliced, corrected, aligned `bed` file
- Sliced, corrected, aligned `reads.fasta` file
- [optional] sliced junctions tab file 

**Inputs for version >= 3**
- `reference.fasta`
- Sliced `reference.gtf`
- Sliced, `bam` file
- [optional] sliced junctions tab file 

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
├── outputs/
│   └── flair_<version>/                         # FLAIR version (e.g., flair_v2.1.1)
│       └── <sample>/                            # Sample identifier (e.g., human_GENCODEv48_WTC11)
│           └── run_<align_id>[_<corr_id>]/      # Alignment settings, Correct settings
│               ├── align/                       # FLAIR align output 
│               │   ├── <sample>.bam
│               │   └── <sample>.bed
│               │
│               ├── correct/                     # FLAIR correct output
│               │   ├── <sample>.corrected.bam
│               │   └── <sample>.corrected.bed
│               │
│               ├── regions/                     # sliced files used for FLAIR collapse runs
│               │   └── chr20_3218000_3250000/
│               │       ├── raw/                 
│               │       │   ├── chr20…-3250000.bam
│               │       │   ├── chr20…-3250000.bed
│               │       │   ├── chr20…-3250000.gtf
│               │       │   ├── chr20…-3250000.fasta
│               │       │   ├── chr20…-3250000.exp5.bed
│               │       │   └── …
│               │       │
│               │       ├── collapse/            # FLAIR collapse per-region, using different options 
│               │       │   └── collapse_<flags>/
│               │       │       ├── *.isoforms.gtf
│               │       │       └── *.isoforms.bed
│               │       │
│               │       ├── collapse_qc/         # FLAIR collapse QC for this region across all runs 
│               │       │   ├── sqanti_summary.tsv
│               │       │   ├── sqanti.png
│               │       │   ├── ted_summary.tsv
│               │       │   ├── ted.png
│               │       │   └── region_metrics.png
│               │       │
│               │       ├── quantify/            # FLAIR quantify per-region, using different options 
│               │       │   └── quantify_<flags>/
│               │       │       ├── isoform_tpm.tsv
│               │       │       └── gene_counts.tsv
│               │       │
│               │       └── quantify_qc/         # FLAIR quantify QC for this region across all runs 
│               │           └── flair_quantify_metrics.png
│               │
│               └── logs/
│                   ├── <align_id>.time.log
│                   ├── <correct_id>.time.log
│                   ├── <collapse_id>.time.log
│                   ├── <quantify_id>.time.log
│                   ├── ted.time.log
│                   ├── sqanti_error.log
│                   └── quantify_qc.log
└── tests/
    └── data/
        └── test_data_set_real/
            ├── sample.gtf
            └── sample_reads.fastq
            └── reference.fastq
            └── [optional] reference TSS/TTS bed file(s)
            └── [optional] experiment TSS/TTS bed file(s)
            └── [optional] reference splice junction bed file
        └── test_data_set_simulated/
            ├── sample.gtf
            └── sample_reads.fastq
            └── reference.fastq
            └── [optional] reference TSS/TTS bed file(s)
            └── [optional] experiment TSS/TTS bed file(s)
            └── [optional] reference splice junction bed file
```


