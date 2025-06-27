# FLAIR Test Suite

A framework for reproducible, version-aware benchmarking and development of FLAIR.

---

## Description 
Compare isoform annotation and quantification across new and exisiting FLAIR releases, sequencing protocols, and genomic regions.

---

## Synopsis
```python
(TBD)
```
---

## Glossary of Terms

**version**
the different versions of FLAIR that are being tested

**flags**
The different options used together at each step of the FLAIR pipeline

**flair_run**
Refers to the entirety of a complete FLAIR run spanning alignment steps to quantification

**QC**
custom and widely used analysis scripts to help us determine the accuracy of each FLAIR run and compre them to each other

**test_regions**
By splitting whole samples into regions - we can run quick tests to assess how different options interact (or how new developments of FLAIR perform relative to past versions).  Regions of interest, derived from the Human GENCODE v48 annotation, are used during development of the test suite. These regions were selected based on qualities such as :
- Region size
- Number of genes
- Total isoforms
- median isoforms per gene
- average gene length
- median of the mean exons per isoform 
- isoform entropy 
- RNA biotypes : such as percentage of RNAs annotated as protein coding vs lncRNA

by selecting on these qualities, we hope to capture the wide range of transcriptome diversity existing within whole samples - while maintaining a balance of speed. 

**test_data_set**
The test data set is what is used for development of the FLAIR test suite. 

`real data`
Attributes: organism, sequencing platform
Data : reads, gene annotation, splice junctions, experimental 5 prime and 3 prime

`simulated data`
Attributes: 

**data_set(s)**
Further data sets can be added to the FLAIR the suite. This includes new organisms, sequencing techniques, simulated reads, and annotations.

**regions**
New regions of interest that can be added to the exisiting regions used in development or exclusively used in the FLAIR test suite

**manifest.py**
file used to set up FLAIR runs including; version, region, flags, QC 

**flair_align_automation.py**
script used to create and run each unique FLAIR align step

**flair_correct_automation.py**
script used to create and run each unique FLAIR correct step

**slice.py**
script used to slice input files by the given regions - to be fed into FLAIR collapse

**flair_collapse_automation.py**
script used to create and run each unique FLAIR collpase step

**sqanti_automation.py**
script used to run sqanti on flair collapse outputs and create the region summary plot 

**ted_automation.py**
script used to analyze transcript end performance on flair collapse outputs and create the region summary plot 

**flair_quantify_automation.py**
script used to create and run each unique FLAIR quantify step

**flair_quantify_qc.py**
(TBD)

---
