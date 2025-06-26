# FLAIR Test Suite Overview

A framework for reproducible, version-aware benchmarking and development of FLAIR.

## 1. Purpose
Compare isoform annotation and quantification across new and exisiting FLAIR releases, sequencing protocols, and genomic regions.

## 2. Getting Started
1. **Data**: drop your FASTA/GTF/BED/BAM files into `tests/data/<dataset_name>/`.  
2. **Config**: edit `manifest.py` to point at your data folder, define regions, and list the FLAIR parameter sets you want to test.

## 3. What It Does
- **Align** reads to the genome  
- **Correct** alignments w/wo splice-junction evidence  
- **Collapse** isoforms into transcript models  
- **Quantify** expression per isoform  

Each stage runs once per parameter set, per region(s)

## 4. Customization & Extension
- **Add more datasets** by creating new `tests/data/<name>/` folders.  
- **Test new FLAIR versions** by appending entries to `align_options`, `correct_options`, etc., each with its own `version`, `env`, and `flags`.  
- **Add new regions** by adding `{ "chrX": ["start:end"] }` to the `regions` list.

---
