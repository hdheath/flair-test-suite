# Overview

## Glossary of Terms

This section defines conventions used across the repo 

(TBD).

---

## Master Data Directory

Where full input data files are stored.

(TBD).

### Expected Input Format

Organize and name experiment folders using this format:

```
<organism>_<reference>_<sample>_<sequencing>_<5pseq>_<3pseq>
```

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

### Adding a New Region

**Before:**

```python
"regions": {
    "chr20": ["3218000:3250000"],
},
```

**After:**

```python
"regions": {
    "chr20": ["3218000:3250000"],
    "chr12": ["4118030:5240000"],
},
```

---

### Adding a New FLAIR Align Option Set

**Before:**

```python
    "align_options": [
        {"--threads": 12},
    ],
```

**After:**

```python
    "align_options": [
        {
            "--threads": 12
        },
        {
            "--threads": 12,
            "--nvrna": True,
        },
    ],
```

### Adding a New FLAIR Correct Option Set

**Before:**

```python
    "correct_options": [
        {
            "--threads": 12,
            "--junc":    os.path.join(DATA_ROOT, "WTC11", "WTC11_all.SJ.out.tab"),
        },
    ],
```

**After:**

```python
    "correct_options": [
        {
            "--threads": 12,
            "--junc":    os.path.join(DATA_ROOT, "WTC11", "WTC11_all.SJ.out.tab"),
        },
        {
            "--threads": 12,
            "--junc":    os.path.join(DATA_ROOT, "WTC11", "WTC11_all.SJ.out.tab"),
            "--ss_window" : 5,
        },
    ],
```

### Adding a New FLAIR Collapse Option Set

**Before:**

```python
"collapse_options": [
    {
        "--version":      2,
        "--support":      1,
        "--check_splice": True,
        "--filter":       "stringent",
        "--threads":      8,
    },
],
```

**After:**

```python
"collapse_options": [
    {
        "version":        2,
        "--support":      1,
        "--check_splice": True,
        "--filter":       "stringent",
        "--threads":      8,
    },
    {
        "version":        3,
        "--support":      3,
        "--filter":       "best-only",
        "--threads":      8,
    },
],
```
