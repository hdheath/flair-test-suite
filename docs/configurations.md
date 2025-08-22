# Guide to Creating Configurations for Test Cases

## **⚠️** Important

* Find config templates in `config/<stages_to_run>_template.toml`.

* All relative paths provided in the config resolve against `data_dir`.

* For targeted region runs:

  * Include the `regionalize` stage in your template.
  * Create a 3‑column TSV file (e.g., `regions.tsv`) with columns: `chr`, `start`, `end`.
  * Save `regions.tsv` in the same directory as your input data.

## Choosing a Template

Select one of the four supported workflows and copy its template from `configs/`:

1. **Align → Correct → Collapse**
2. **Align → Regionalize → Correct → Collapse**
3. **Align → Transcriptome** (requires **FLAIR ≥ 3.0**)
4. **Align → Regionalize → Transcriptome** (requires **FLAIR ≥ 3.0**)

## Editing and Saving Templates

Once you have selected and copied a template, you can edit it by : 

1. Set `run_id` under `[run]` to name your output folder.
2. Edit the six fields under `[run]`: `version`, `conda_env`, `work_dir`, `data_dir`, `reads_file`, `genome_fa`.
3. Add any CLI flags under each `[run.stages.flags]`; leave blank to use FLAIR defaults.
  **⚠️** Do not use -o flag

4. In addition to providing input data, users should fill a description which provides reasoning as to why a test case is created. This description should be commented in the configuration file, and would include things like : 

    | Field              | Description                                                                                                           |
    |--------------------|-----------------------------------------------------------------------------------------------------------------------|
    | `organism`         | Organism from which reads are derived                                                          |
    | `sample(s)`         | Sample(s) from which reads are derived
    | `assembly`         | Version of reference genome assembly                                                                                  |
    | `platform`         | Sequencing platform (e.g. ONT, PacBio, Simulated)                                                                     |
    | `region selection` | Reasons for selecting targeted region(s) for testing:<br><ul><li>Region size</li><li>Number of genes</li><li>Total isoforms</li><li>Median isoforms per gene</li><li>Average gene length</li><li>Median of mean exons per isoform</li><li>Isoform entropy</li><li>RNA biotypes (e.g., % protein coding vs. lncRNA)</li></ul> |

5. Save your edited TOML file in the `configs/` directory.

## Example of an Edited Configuration

```toml
# ------------------------------------------------------------------------------------------------------------
# Align -> Correct -> Regionalize -> Collapse Configuration
# ------------------------------------------------------------------------------------------------------------
# -------------------------------------------- Test Case Description -----------------------------------------
# 
# - Organism: Human
# - Sample(s) : WTC11
# - Assembly: GR38
# - Platform: PacBio

# -------------------------------------------- Region Selection 
# chr5:14000-17000 : single gene, w lots of exons per isoform
# chr12:3500-160000 : lncRNA

# -------------------------------------------- Summary
# Testing default collapse options on lncRNA and single gene region

# ------------------------------------------------------------------------------------------------------------

run_id = "WTC11"

[run]
version    = "3.0.0"                # FLAIR version
conda_env  = "flair"                # conda env name
work_dir   = "./outputs"            # output directory
data_dir   = "./tests/data"         # input data directory
reads_file = "WTC11_reads.fasta"    # long-read FASTA
genome_fa  = "gr38.gtf"             # reference genome FASTA

# ------------------------------------------------------------------------------------------------------------
# -------- ALIGN ------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

[[run.stages]]
name = "align"

[run.stages.flags]
nvrna   = true  # on/off flag
quality = 0

# ------------------------------------------------------------------------------------------------------------
# -------- CORRECT ----------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

[[run.stages]]
name     = "correct"
requires = ["align"]

[run.stages.flags]
gtf          = "hs_GENCODE38.basic_annotation.gtf"
junction_bed = "WTC11_all.SJ.out.tab"

# ------------------------------------------------------------------------------------------------------------
# -------- REGIONALIZE ----------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

[[run.stages]]
name     = "regionalize"
requires = ["correct"]

[run.stages.flags]
gtf                                 = "hs_GENCODE38.basic_annotation.gtf"
regions_tsv                         = "regions_of_interest.tsv"
junctions                           = "WTC11.tab"                 # optional
experiment_5_prime_regions_bed_file = ""                          # optional
experiment_3_prime_regions_bed_file = ""                          # optional
reference_5_prime_regions_bed_file  = ""                          # optional
reference_3_prime_regions_bed_file  = ""                          # optional

# ------------------------------------------------------------------------------------------------------------
# -------- COLLAPSE ----------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

[[run.stages]]
name     = "collapse"
requires = ["regionalize"]

[run.stages.flags]
                  # Note : No flags, so default settings for collapse used
```

## Combine stage manifest

The `combine` stage merges multiple transcriptomes using a manifest TSV.
Supply an existing manifest file via the `manifest` flag and the collapse
output from the current run will be appended to the end of that file. If no
manifest is provided, only the collapse output is used.

```toml
[[run.stages]]
name = "combine"
requires = ["collapse"]

[run.stages.flags]
manifest = "my_manifest.tsv"
```

Each manifest row has five tab-separated columns:

``sample`` `type` `bed` `fasta` `read_map`

The `sample` field must be a numeric identifier. When a manifest is
provided, collapse outputs from the current run are appended using the next
available sample number (maximum existing `sample` + 1). Paths are resolved
relative to `data_dir`. Leave `fasta` and `read_map` blank if those files are
unavailable.

For more details, see the [FLAIR Test Suite Overview](./overview.md).

## Quantify stage manifest

The `quantify` stage computes isoform and gene expression using a reads
manifest. It can follow `collapse`, `transcriptome`, or `combine`. A manifest
file is **required** to list all FASTQ files used for expression quantification.

```toml
[[run.stages]]
name = "quantify"
requires = ["collapse"]   # or ["transcriptome"] or ["combine"]

[run.stages.flags]
manifest = "reads_manifest.tsv"
```

Each manifest row has four tab-separated columns:

``sample`` `condition` `batch` `reads.fq`

Entries for the current run's `reads_file` (or `reads_files`) are appended
automatically. Sample names must be unique and paths should point to existing
FASTQ files, preferably using absolute paths.
