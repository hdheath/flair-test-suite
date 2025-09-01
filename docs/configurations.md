# Guide to Creating Configurations for Test Cases

## Table of Contents

- Quick Start
- Combined TSV Format
- Stage Ordering Rules
- Behavior & Effects (Overview)
- Config Cookbook

## Quick Start

1) Copy and edit the configuration template :

```bash
cp config/templates/flair_test_suite_config.tsv name/your/config.tsv
```

- Fill the header on the first `align` row: `test_set_id`, `version`, `flair_env`, `data_dir`, `reads_file`, `genome_fa` (as needed), and shared inputs like `gtf`.
- Add stage rows in order; put CLI flags in the `flags` column (comma-separated tokens).

2) Run the suite:

```bash
flair-test-suite flair_test_suite_config.tsv
```

3) Verify outputs:

- Results: `./outputs/<test_set_id>/`
- Log: `./outputs/<test_set_id>/run_summary.log`

## Combined TSV Format

- Command: `flair-test-suite <flair_test_suite_config.tsv>`
- One combined file per test set; multiple cases are supported. A new case starts at an `align` row.
- Header columns (case-insensitive):
  - `test_set_id`, `version`, `flair_env`, `data_dir`, `reads_file`, `genome_fa`, `gtf`, optional `junctions`, `sqanti_env`, and required `stage`, `flags`.
- Flags: enter the exact CLI flags you would pass to FLAIR, but comma-separated within the cell. Both `--opt=value` and `--opt value` forms are accepted. Quote values that contain commas.
- The harness will automatically enforce required options: `-r/--reads`, `-g/--genome`, `-q/--bed`, `-b/--bam`, `-o/--out`.
- Region runs: use stage name `region_test` and pass `--regions-tsv=regions.tsv` in the flags column. The regions TSV has 3 columns: `chr  start  end`.

### Short Example (single case)

`flair_test_suite_config.tsv`

```text
key	value
test_set_id	WTC11_demo
version	3.0.0
flair_env	flair
data_dir	/data/WTC11
reads_file	reads.fastq.gz
genome_fa	GRCh38.fa
gtf	gencode.v48.gtf
experiment_5_prime_regions_bed_file	exp_TSS.bed
experiment_3_prime_regions_bed_file	exp_TES.bed
reference_5_prime_regions_bed_file	ref_TSS.bed
reference_3_prime_regions_bed_file	ref_TES.bed

stage	flags
align	--nvrna, --threads 8
region_test	--regions-tsv regions.tsv
correct	--gtf
collapse	
```


### Multiple Cases (cache reuse)

You can define multiple cases in the same configuration tsv. Reused stages with identical inputs/flags are skipped and their outputs are reused (via the stage signature under `outputs/<test_set_id>/<stage>/<signature>/`).

```text
key	value
test_set_id	WTC11_demo
version	3.0.0
flair_env	flair
data_dir	/data/WTC11
reads_file	reads.fastq.gz
genome_fa	GRCh38.fa
gtf	gencode.v48.gtf

stage	flags
# Case 1 — collapse path
align	--nvrna, --threads 8
correct	
collapse	

# Case 2 — transcriptome path (reuses cached align)
align	--nvrna, --threads 8
transcriptome	
```

### Reads input forms

- Single file: `reads_file = "reads.fastq.gz"`
- List multiple files as list-like string with commas: `reads_file = "r1.fq.gz,r2.fq.gz"`

## Stage Ordering Rules

- `region_test` must follow `align`.
- `correct` must follow `align` (optionally after `region_test`) and cannot follow `collapse`, `transcriptome`, `combine`, or `quantify`.
- `collapse` must follow `correct`.
- `transcriptome` must follow `align` (optionally after `region_test`) and cannot follow `correct`.
- `combine` must follow `collapse` or `transcriptome` unless `--manifest` is provided.
- `quantify` must follow `combine`, `transcriptome`, or `collapse`.

## Behavior & Effects (Overview)

For output locations, re‑run and caching rules, and region‑test behavior, see the Overview:

- [Overview: Output Layout](overview.md#output-layout)
- [Overview: Re‑run & Caching](overview.md#re-run--caching)
- [Overview: Regionalized Effects](overview.md#regionalized-effects)

## Config Cookbook

- Targeted region runs: add `region_test` after `align`; pass `--regions-tsv=regions.tsv` in the flags column.
- Multiple reads: use a comma‑string in `reads_file`.
- TED for collapse/transcriptome runs automatically; provide peak BEDs under the header to enable precision/recall metrics; otherwise those metrics are `None`.
- SQANTI runs automatically when a suitable conda env is available. Set `sqanti_env` in the header; CPUs are fixed to 4.
- Combine with optional manifest: include a `combine` row and pass `--manifest=manifest.tsv`.
- Quantify with reads manifest: include a `quantify` row and pass `--manifest=reads_manifest.tsv`.
