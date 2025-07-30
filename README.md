# FLAIR Test Suite 🧪

*A config‑driven regression and QC harness for*
[FLAIR](https://github.com/BrooksLabUCSC/flair) **≥ 2.0** long‑read transcriptome analysis.

> **Highlights**
> • Stage‑level caching (skip work when nothing changed)
> • Automated QC checkpoints after every stage (PNG + TSV)

---

## Table of Contents

1. [Design](#design)
2. [Dataset Types](#dataset-types)
3. [Region / Scope](#regionscope)
4. [Installation & Setup](#installation--setup-single‑step-conda--pip)
5. [Configuration Guide](#configuration-guide)
6. [Running the Test Suite](#running-the-test-suite)
7. [Pipeline Stages & QC](#pipeline-stages--qc)
8. [Output Layout](#output-layout)
9. [Contributing & License](#contributing--license)

---

## Design

The suite is organised around **end‑to‑end test cases**.
Each test case runs a full FLAIR workflow defined by a single TOML file.
After every stage, QC scripts evaluate outputs. 

*Stage signature* = `tool_version | flags | input_hashes`
If the signature hasn’t changed, the stage is skipped and outputs are reused.

---

## Dataset Types

| Type          | Purpose                                                                                      |
| ------------- | -------------------------------------------------------------------------------------------- |
| **Simulated** | Truth isoforms known → rigorous correctness checks                                           |
| **Real data** | ONT / PacBio (cDNA, direct RNA, high/low depth) → tests platform and library‑specific issues |

*(Dataset attributes: name, description, data\_source, organism, assembly, platform, read\_fastq, gene\_set, optional TSS/TES/junction evidence …)*

---

## Region/Scope

| Scope                   | Use‑case                                    |
| ----------------------- | ------------------------------------------- |
| **Targeted regions**    | Fast, focused edge‑case tests (single loci) |
| **Whole‑transcriptome** | Scalability / performance validation        |

*(Region attributes: chr, start, end; considerations: gene count, isoform entropy, etc.)*

---

## Installation & Setup

> Requires Conda or Mamba.

```bash
git clone https://github.com/hdheath/flair-test-suite.git
cd flair-test-suite

# create env, install FLAIR + samtools + this package (+ all Python deps)
conda env create -f flair-test-suite.yaml
conda activate flair-test-suite
```

### Verify

```bash
which flair-test-suite        # should print path inside env
```

---

## Configuration Guide

*Important*

* Users must already have downloaded the FLAIR conda env they would like to use.
* All data used for a run must be located in the same directory 

* Users choose to run from 4 full FLAIR workflows
    1. Align -> Correct -> Collapse
    2. Align -> Correct -> Slice -> Collapse
    3. Align -> Transcriptome
    4. Align -> Slice -> Transcriptome

*Key Points*
* Each configuration file must be created by the user;
    * Find run templates in **`configs/<stages_to_run>_template.toml`**.
    * All relative paths provided in the config resolve against `data_dir`.
    * Customize FLAIR stage flags by adding them under each `[run.stages.flags]`. (otherwise defaults used)

* For Targeted Region Runs;
    * Make sure your config includes the **slice** stage
    * create a 3 column tsv (eg <regions.tsv>) in the following format : (**chr**\t**start**\t**end**)
    * save <regions_tsv> under the directory where your input data is 

### Example Configuration File
```toml
# ------------------------------------------------------------------------------------------------------------
# Align -> Correct -> Slice -> Collapse Configuration
# ------------------------------------------------------------------------------------------------------------

run_id = "WTC11"

[run]
version    = "3.0.0"                # version of FLAIR
conda_env  = "flair"                # name of user installed conda FLAIR env
work_dir   = "./outputs"            # where output will be saved
data_dir   = "./tests/data"         # where data lives for run
reads_file = "WTC11_reads.fasta"    # reads fasta
genome_fa  = "gr38.gtf"             # reference genome fasta

# ------------------------------------------------------------------------------------------------------------
# -------- ALIGN ------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

[[run.stages]]
name = "align"

[run.stages.flags]
nvrna   = true  # set on/off flags to true
quality = 0

# ------------------------------------------------------------------------------------------------------------
# -------- CORRECT ----------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
[[run.stages]]
name     = "correct"
requires = ["align"]

[run.stages.flags]
gtf              = "hs_GENCODE38.basic_annotation.gtf"  
junction_bed     = "WTC11_all.SJ.out.tab"

# ------------------------------------------------------------------------------------------------------------
# -------- SLICE ----------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
[[run.stages]]
name = "slice"
requires = ["correct"]

[run.stages.flags]
gtf = "hs_GENCODE38.basic_annotation.gtf" 	
regions_tsv = "regions_of_interest.tsv"
junctions = "WTC11.tab"   			              # optional
experiment_5_prime_regions_bed_file = ""  # optional
experiment_3_prime_regions_bed_file = ""  # optional
reference_5_prime_regions_bed_file =  ""  # optional
reference_3_prime_regions_bed_file =  ""  # optional

# ------------------------------------------------------------------------------------------------------------
# -------- COLLAPSE ----------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
[[run.stages]]
name = "collapse"
requires = ["slice"]

[run.stages.flags]

```

---

## Running the Test Suite

You can run individual test‑cases or batch multiple configurations.

### Single config

```bash
flair-test-suite path/to/your_config.toml
```

This will execute the full pipeline defined in `your_config.toml` and place results under `<work_dir>/<run_id>`.

### Batch runs from a list

If you have many configs, create a text file listing each TOML path, one per line. For example, `runs.txt`:

```
configs/demo.toml
configs/sample1.toml
configs/sample2.toml
```

Then launch them in sequence:

```bash
flair-test-suite path/to/runs.txt
```

---


## Pipeline Stages & QC

| Stage           | Primary output(s)     | QC metrics / plots                                          |
| --------------- | --------------------- | ----------------------------------------------------------- |
| `align`         | BAM + BED             | MAPQ, read identity/length, unique junctions, splice motifs |
| `correct`       | corrected BED         | reads removed %, unique junctions, splice‑motifs    |
| `slice`         | region BAM/BED/FA/GTF | feature counts from gtf  (such as number of genes)                             |
| `collapse`      | isoforms BED/GTF      | *QC TBD*                                                    |
| `transcriptome` | isoforms BED/GTF      | *QC TBD*                                                    |

PNG plots and TSV metrics are saved next to each stage’s outputs.

---

## Output Layout

```plaintext
<work_dir>/
└── <run_id>/
    ├── align/<sig>/
    │   ├── <run_id>_flair.bam
    │   ├── <run_id>_flair.bed
    │   └── align_qc.{tsv,png}
    ├── correct/<sig>/
    ├── slice/<sig>/
    ├── collapse/<sig>/
    └── transcriptome/<sig>/
```

---

## Contributing & License

* Bug reports and PRs welcome — see **`CONTRIBUTING.md`**.
* Code of Conduct: **`CODE_OF_CONDUCT.md`**.
* © 2025 **Harrison Heath / Brooks Lab** – released under the **MIT License** (`LICENSE`).

Happy testing 🚀
