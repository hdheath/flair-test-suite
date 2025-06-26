# Overview

## Glossary of Terms

This section defines conventions used across the repo 

(TBD).

---

## Master Data Directory

All full input datasets should be stored under:

(TBD).

---

### Expected Directory Name Format

Each experiment should follow a naming convention that includes organism, reference annotation, sample name, sequencing platform, and optional 5' or 3' protocol descriptors.

```
<organism>_<reference>_<sample>_<sequencing>_<5pseq>_<3pseq>
```
---

### Example Data Directory Layout

```plaintext
flair-test-suite/
└── tests/
    └── data/
        └── human_GENCODEv48_WTC11_PacBio_CAGE_QuantSeq/
            ├── sample.gtf
            ├── sample_reads.fastq
            ├── reference.fastq
            ├── reference TSS/TTS BED file(s)
            ├── experiment TSS/TTS BED file(s)
            └── reference splice junction BED file
        └── human_GENCODEv48_WTC11_Nanopore/
            ├── sample.gtf
            ├── sample_reads.fastq
            ├── reference.fastq
            ├── reference TSS/TTS BED file(s)
            └── reference splice junction BED file
```

---

## Editing the Manifest

The manifest file defines what regions and pipeline settings to use. You can extend its scope by adding new entries to any of the following sections.

---

### Adding a New Region

Each chromosome key maps to a list of coordinate ranges in the format "chr:start-end". Each cooridinate range will be processed independently.

To add a new genomic region for processing :
- Locate the `regions` section in your manifest file
- Add a new dictionary to the existing coordinate range list of dictionaries with your desired range. (ex : `{chr1:140000-150000}`)

---

### Adding a New FLAIR Option Set (Align, Correct, Collapse, Quantify)

The manifest includes separate lists for different stages of the FLAIR pipeline: `align_options`, `correct_options`,`collapse_options`, `quantify_options`. Each list contains dictionaries, where each dictionary specifies a distinct set of flags for that step.

To add a new configuration to any of these stages:
- Append a new dictionary with the desired flags to the appropriate list.
- Each configuration will trigger an independent run of the relevant FLAIR stage.

---