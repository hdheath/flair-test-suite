# FLAIR Test Suite Overview

The FLAIR Test Suite is organized to cover the pipeline’s functionality across multiple dimensions. Test cases are grouped in ways that make it easy to target specific aspects of the pipeline or data contexts. The primary groupings are:

---

### Design

The FLAIR Test Suite is organized around **end-to-end test cases**, where each **test-case** is a complete FLAIR workflow from raw reads through FLAIR quantify with self-contained stage options at each step. After each stage of a test-case, a set of **QC Checkpoints** (sub-tests) validate the intermediate outputs.

---

### Dataset Type

Test cases are defined by the nature of the dataset:

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

### Region/Scope

Test cases are also defined by region. A region is a genomic coordinate where we expect FLAIR to analyze transcripts. Test cases choose one of these region scopes:

- **Targeted Region Tests**  
  Run on a limited locus or small gene set (e.g., reads mapping to chr21 or a single gene). Cover challenging regions (e.g., high gene density, pseudogenes, repetitive sequences). Helps ensure corner cases are regularly checked.

- **Whole-Transcriptome Tests**  
  Run on genome-wide data (e.g., whole human transcriptome) to ensure the pipeline scales to full dataset sizes and complexities.

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

### FLAIR Version

In each test case you specify the **FLAIR version** to run (for example, `v1.5` or `v2.0-dev`). 


---

