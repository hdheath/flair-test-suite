## FLAIR Test Suite Overview

The FLAIR Test Suite is organized to cover the pipeline’s functionality across multiple dimensions. Test cases are grouped in ways that make it easy to target specific aspects of the pipeline or data contexts. The primary groupings are:

---

### By Workflow Mode

The FLAIR Test Suite is organized around **end-to-end test cases** and **stage-specific test cases**.

#### 1. End-To-End Test Cases  
Run the complete FLAIR workflow from raw reads through QC (and optional differential analysis) with self-contained stage options at each step. After each stage of an end-to-end test case, a set of **QC Checkpoints** (sub-tests) validate the intermediate outputs.

#### 2. Stage Test Cases  
Run a specific stage of the FLAIR pipeline (e.g., align, correct, collapse, …), from a well-defined input state. Stage tests are categorized by the module they exercise. Grouping by stage allows testers to run a focused subset of tests. For instance, during development of an improved collapse algorithm, one can run just the Collapse Stage Tests to verify its behavior without running the full suite.

---

### By Dataset Type

Test cases are also grouped by the nature of the dataset:

- **Simulated Data Tests**  
  Use artificial datasets where the “true” isoforms are known in advance (e.g., simulated reads from a known transcript set, or spike-in controls). These help validate correctness without biological ambiguity.

- **Real Experimental Data Tests**  
  Use real sequencing datasets (e.g., human nanopore cDNA, direct RNA, PacBio Iso-Seq). While ground truth is not fully known, these tests come with expected biological behaviors (e.g., known mutation effects or tissue-specific splicing patterns).

  - **Platform/Library Variants**  
    Group cases by sequencing platform or library prep (e.g., Oxford Nanopore vs PacBio; cDNA vs direct RNA). Also test different read lengths or depths (e.g., high-depth vs low-depth) to ensure FLAIR performs robustly.


This grouping ensures that any changes affecting a specific data type (e.g., poly(A) tail handling) can be checked in isolation, and helps identify data-specific issues.

#### *Dataset* attributes
- *name*  
- *description*  
- *data_source* – real or simulated  
- *organism*  
- *assembly*  
- *platform* – ONT or PB  
- *read_fastq* – one or more FASTQ files  
- *gene_set_name*  
- *gene_set_annot_gtf* – GTF annotation file  
- *gene_set_meta_tsv* – TSV with metadata (extracted from GTF)  
- *splice_junctions_files* – splice junction files (STAR, BED)  
- *tss_evidence* – TSS evidence files (CAGE, etc)  
- *tes_evidence* – TES evidence files (QuantSeq, etc)  
- *ground_truth* – for simulated data  

---

### By Region/Scope

A region is a genomic coordinate where we expect FLAIR to analyze transcripts. Test cases choose one of these region scopes:

- **Targeted Region Tests**  
  Run on a limited locus or small gene set (e.g., reads mapping to chr21 or a single gene). These quick smoke tests are useful for debugging specific issues.

- **Whole-Transcriptome Tests**  
  Run on genome-wide data (e.g., whole human transcriptome) to ensure the pipeline scales to full dataset sizes and complexities.

- **Regional Edge-Case Tests**  
  Cover challenging regions (e.g., high gene density, pseudogenes, repetitive sequences). Helps ensure corner cases are regularly checked.

Grouping by region helps select quick tests for fast feedback or full-scale validation to reveal scale-dependent issues (e.g., memory leaks).

#### *Region* attributes
- *chr*  
- *start*  
- *end*  

#### *Region* considerations for selection
- Region size  
- Number of genes  
- Total isoforms  
- Median isoforms per gene  
- Average gene length  
- Median of mean exons per isoform  
- Isoform entropy  
- RNA biotypes (e.g., % protein coding vs lncRNA)  

---

### By FLAIR Version

In each test case you specify the **FLAIR version** to run (for example, `v1.5` or `v2.0-dev`). 


---

## Summary of Test Suite Structure

In practice, each test case is uniquely identified and tagged along four dimensions:

| Dimension   | Example Tags                             |
|-------------|------------------------------------------|
| **Workflow**| `E2E`, `AlignOnly`, `CollapseOnly`       |
| **Dataset** | `Simulated`, `PacBio HiFi`               |
| **Region**  | `chr21:42000-62000,chr13:14300-15000`, `WholeTranscriptome`, `chr13` |
| **Version** | `v1.5`, `v2.0-dev`, `v2.0-newFeatures`   |
 

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

