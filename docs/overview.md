# FLAIR Test Suite

A framework for reproducible, version-aware benchmarking and development of FLAIR.

---

## Description 
Compare isoform annotation and quantification across new and exisiting FLAIR releases, sequencing protocols, and genomic regions.

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
Regions of interest, derived from the Human GENCODE v48 annotation, and used during development of the test suite. These regions were selected based on qualities such as :
Region size
Number of genes
Total isoforms
median isoforms per gene
average gene length
median of the mean exons per isoform 
isoform entropy 
RNA biotypes : such as percentage of RNAs annotated as protein coding vs lncRNA

**test_data_set**
The test data set is what is used for development of the FLAIR test suite. This includes :

`real data`
Input data : 
Organism : Human
Sequencing : PacBio
Sample : WTC11
Annotation : GENCODE v48
Splice Junctions : 

Data for QC : 
Experimental 5 prime data : CAGE
Experimental 3 prime data : QuantSeq

`simulated data`
Input data : 
Sequencing : Simulated long reads 

Data for QC : 
5 prime and 3 prime data : simulated read ground truth 

**data_set(s)**
Further data sets can be added to perform QC of FLAIR. This includes new organisms, sequencing techniques, simulated reads, and annotations.

**regions**
New regions of interest that can be added to the exisiting regions used in development or exclusively used in the FLAIR test suite

---
