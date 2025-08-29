# Guide to Creating Configurations for Test Cases

## Table of Contents

 - [Quick Start: Create & Run](#quick-start-create--run)
 - [Authoring a Test Case](#authoring-a-test-case)
 - [Templates Overview](#templates-overview)
 - [TOML Schema & Keys](#toml-schema--keys)
 - [Data Input Example](#data-input-example)
 - [Flags Syntax](#flags-syntax)
   - [Flags → CLI mapping](#flags--cli-mapping)
   - [Reads input forms](#reads-input-forms)
 - [CLI Reference (TSV-only)](#cli-reference-tsv-only)
 - [Create a New Run (TSV-based)](#create-a-new-run-tsv-based)
   - [Common Pitfalls & Fixes](#common-pitfalls--fixes)
   - [Build the TSV (examples)](#build-the-tsv-examples)
   - [Optional test-case notes](#optional-test-case-notes)
 - [Stage Ordering Rules](#stage-ordering-rules)
 - [Output Layout & Caching](#output-layout--caching)
 - [Config Cookbook](#config-cookbook)
 - [Tips](#tips)
  - [Regionalized Effects](#regionalized-effects)

## Quick Start: Create & Run

1) Copy and edit base inputs:

```bash
cp config/templates/inputs.toml config/my_inputs.toml
```

- Set `test_set_id` and fill `[run]` fields (`version`, `conda_env`, `data_dir`, `reads_file`, and any shared paths like `gtf`).

2) Pick a workflow and copy stage templates you need (edit flags under `[run.stages.flags]` as desired):

- Collapse path: `align.toml` → `correct.toml` → `collapse.toml`
- Transcriptome path (FLAIR ≥ 3.0): `align.toml` → `transcriptome.toml`
- Targeted runs: insert `regionalize.toml` after `align.toml`

3) Build a TSV manifest (one test case per line):

```text
my_inputs.toml	config/templates/align.toml	config/templates/correct.toml	config/templates/collapse.toml
```

4) Run the suite:

```bash
flair-test-suite path/to/cases.tsv
```

5) Verify outputs:

- Results: `./outputs/<test_set_id>/`
- Log: `./outputs/<test_set_id>/run_summary.log`



## Important

- The test suite runs from a TSV configuration file (template: `config/templates/test_case_cli_template.tsv`).
  - Each row is a test case.
  - A test case lists the stages to run in order.
  - Each stage has its own TOML configuration.
- Stage templates live under `config/templates/<stage_name>/`.
- Reuse stage configs across multiple test cases.
  - Example: one `align.toml` can feed different downstream workflows.
  - The suite reuses cached outputs from upstream stages automatically.
- Relative input paths resolve against `run.data_dir`; template paths in the TSV resolve against the TSV’s location.
- To run test cases on a user-subset genomic region : 
  - Include a `regionalize` stage.
  - Provide a 3‑column TSV (`chr`, `start`, `end`), e.g., `config/templates/region_template.tsv` specifying the genomic range(s) you would like your data to be subset to (one range per line)

## Authoring a Test Case

This guided flow mirrors how you actually build a runnable case: base inputs → stage configs → TSV manifest → run.

Step 1 — Base inputs (required)

1) Copy the base template and edit values:

```bash
cp config/templates/inputs.toml my_inputs.toml
```

2) In `my_inputs.toml`, set:
- `test_set_id`
- `[run]`: `version`, `conda_env`, `data_dir`, `reads_file`, `genome_fa` (as applicable), plus shared inputs like `gtf`, `regions_tsv`.
- Optional: provide TED peak BEDs under `[run]` to enable TED precision/recall metrics (TED runs automatically for collapse/transcriptome; without peaks these metrics are `None`).

Reference: see “TOML Schema & Keys” and “Data Input Example” below.

Step 2 — Stage templates (compose workflow)

1) Choose a path:
- Collapse: align → correct → collapse
- Transcriptome (FLAIR ≥ 3.0): align → transcriptome
- Targeted: insert regionalize after align

2) Edit flags in each stage file under `[run.stages.flags]` as needed:
- `config/templates/align.toml`
- `config/templates/correct.toml`
- `config/templates/regionalize.toml` (optional)
- `config/templates/collapse.toml` or `config/templates/transcriptome.toml`
- Optional: `config/templates/combine.toml`, `config/templates/quantify.toml`

3) Optional QC:
- TED: runs automatically for `collapse` and `transcriptome`. Peaks in `[run]` are optional; when absent, precision/recall are `None`.
- SQANTI: runs automatically when a conda env is available. Prefer setting `[run].sqanti_env`; if not set, the harness uses `sqanti` when it exists. Use `[qc.<stage>.SQANTI]` only to override CPUs per-stage.

Step 3 — Build the TSV and run

1) Create a TSV where each line is a case: base inputs first, then stage TOMLs left→right in execution order. Paths are resolved relative to the TSV file.

Example (collapse workflow):

```text
my_inputs.toml	config/templates/align.toml	config/templates/correct.toml	config/templates/collapse.toml
```

2) Run:

```bash
flair-test-suite path/to/cases.tsv
```

3) Check outputs/logs:
- Outputs: `./outputs/<test_set_id>/`
- Log: `./outputs/<test_set_id>/run_summary.log`

Notes:
- TSV order defines upstreams and is validated (“Stage Ordering Rules”).
- Outputs always write to `./outputs/`.

## CLI Reference (TSV‑only)

- Command: `flair-test-suite <cases.tsv>`
- Format (tab‑separated, paths resolve relative to the TSV file):

```text
base_config.toml	stage1.toml	stage2.toml	...
```

Notes:
- One test case per line; stage order is left→right.
- Outputs are always written under `./outputs/<test_set_id>/`.

## TOML Schema & Keys

- `test_set_id` (string): required identifier for the run. Becomes the folder under `./outputs/` and the output filename prefix.
- `[run]` (table): required run‑level inputs shared by stages.
  - `version` (string): FLAIR version, e.g., `"3.0.0"`.
  - `conda_env` (string): conda environment containing FLAIR (and optional SQANTI).
  - `data_dir` (string): base directory for input data. All relative file paths resolve against this directory.
  - `reads_file` (string | list[string]): single path, comma‑separated string (`"r1.fq,r2.fq"`), or list. Escape commas within a path using `\,`.
- `genome_fa` (string): reference genome FASTA.
  - Optional shared inputs: `gtf`, `junctions`, `regions_tsv`.
  - Optional tool settings:
    - `sqanti_env` (string): preferred conda env name for SQANTI (e.g., `"sqanti"`). When set or when an env named `sqanti` exists, SQANTI runs automatically for collapse/transcriptome.
    - `sqanti_cpus` (int): CPUs for SQANTI (optional; default 4).
  - Optional TED peak BEDs (enable precision/recall metrics for TED; TED runs automatically):
    - `experiment_5_prime_regions_bed_file`, `experiment_3_prime_regions_bed_file`
    - `reference_5_prime_regions_bed_file`, `reference_3_prime_regions_bed_file`
- `[[run.stages]]` (array): ordered stage entries (one per template you include in the TSV line).
  - `name` (string): `align`, `regionalize`, `correct`, `collapse`, `transcriptome`, `combine`, or `quantify`.
  - `flags` (table, optional): stage‑specific CLI flags.
- `[qc]` (table): optional QC blocks merged from any stage fragment you include.
  - TED: no block required; TED runs automatically for `collapse` and `transcriptome` with window fixed to 50.
  - SQANTI: runs automatically if a suitable conda env is available. Prefer `[run].sqanti_env` to select the environment. `[qc.<stage>.SQANTI]` is optional for per-stage CPU overrides.

## Data Input Example

Begin by creating the input configuration file from the template at `config/templates/inputs.toml`. These run‑level inputs are used throughout the entire test case.

```toml
test_set_id = "WTC11_demo"

[run]
version = "3.0.0"
conda_env = "flair"
data_dir = "/data/WTC11"
reads_file = "reads.fastq.gz"
genome_fa = "GRCh38.fa"
gtf = "gencode.gtf"
# Optional for regionalize
regions_tsv = "regions.tsv"

# Optional: enable TED QC by providing these paths
experiment_5_prime_regions_bed_file = "exp_TSS.bed"
experiment_3_prime_regions_bed_file = "exp_TES.bed"
reference_5_prime_regions_bed_file  = "ref_TSS.bed"
reference_3_prime_regions_bed_file  = "ref_TES.bed"
```

## Templates Overview

Templates live in `config/templates/` and are referenced in the TSV. Pick what you need:

- `align.toml`: aligns reads to `run.genome_fa`; outputs `<id>_flair.{bam,bed}`.
- `regionalize.toml` (optional): subsets to regions in `run.regions_tsv`; writes per‑region BAM/BED/GTF and `qc/regionalize/region_details.tsv`.
- `correct.toml`: consumes align/regionalize BED; writes `<id>_all_corrected.bed` or `{tag}_all_corrected.bed`.
- `collapse.toml` or `transcriptome.toml`: produces `{*.isoforms.{bed,gtf}}` (per‑region when regionalized). TED QC runs automatically; provide peak BEDs under `[run]` to enable precision/recall metrics.
- `combine.toml` (optional): merges isoforms; auto‑discovers upstream isoforms and appends to manifest.
- `quantify.toml` (optional): quantifies from upstream isoforms; requires a reads manifest TSV.

## Flags Syntax 

Stage-specific options live under `[run.stages.flags]` in each stage’s TOML.

- Presence-only flags: set to true/empty/null to emit the flag
  - `nvrna = true` → `--nvrna`
- Key/value flags:
  - Numbers: `threads = 8` → `--threads 8`
  - File paths: `bed = "regions.bed"` → resolves relative to `run.data_dir` (or used as-is if absolute) and is included in the signature
- Do not set output names (`-o`/`--out`); the harness manages outputs and directories.

Note: Bare keys (e.g., writing just `nvrna` without `= value`) are not valid TOML and are not supported.

### Flags → CLI mapping

- Presence flags: set value to `true`, empty string `""`, or `null` → emit just the flag
  - `nvrna = true` → `--nvrna`
- Number values: `threads = 8` → `--threads 8` (or `-t 8` if you use key `t`)
- String values:
  - If the value resolves to an existing file under `run.data_dir`, it is emitted as a file path and included in the stage’s signature (for caching).
  - Otherwise the string is passed through literally as the flag value.
- File flags are hashed for caching: changing a referenced file or flag value changes the signature and creates a new stage directory.

Examples (inside `[run.stages.flags]`):

```toml
# Presence only
nvrna = true            # -> --nvrna

# Short vs long; both work
t = 8                   # -> -t 8
threads = 8             # -> --threads 8

# File paths (resolve against run.data_dir unless absolute)
bed = "regions.bed"     # -> --bed /abs/path/to/regions.bed
index = "/abs/index"   # -> --index /abs/index

# Literal strings (no file at value)
mode = "speed"          # -> --mode speed
```

### Reads input forms

- Single file: `reads_file = "reads.fastq.gz"`
- Multiple files as list: `reads_file = ["r1.fq.gz", "r2.fq.gz"]`
- Multiple files as one string: `reads_file = "r1.fq.gz,r2.fq.gz"`
  - Escaped commas inside a path: `"dir\\,with\\,comma/sample.fq"`

## Stage Templates (Quick Reference)

- align.toml
  - Needs: `run.genome_fa`, `run.reads_file`
  - Outputs: `<test_set_id>_flair.{bam,bed}`; QC under `qc/`
- regionalize.toml (optional)
  - Needs: `run.gtf`, `run.regions_tsv` (3‑col TSV)
  - Outputs: per‑region BAM/BED/GTF; index at `qc/regionalize/region_details.tsv`
- correct.toml
  - Consumes: align/regionalize BED
  - Outputs: `<id>_all_corrected.bed` or `{tag}_all_corrected.bed`; per‑region QC
- collapse.toml or transcriptome.toml
  - Produces isoforms `{*.isoforms.{bed,gtf}}` (per‑region when regionalized)
  - TED QC runs automatically; add peak BEDs under `[run]` for precision/recall
  - SQANTI runs automatically when `[run].sqanti_env` is set (or `sqanti` env exists)
- combine.toml (optional)
  - Auto‑discovers upstream isoforms and appends to manifest; outputs merged BED/FA/Counts
- quantify.toml (optional)
  - Requires a reads manifest TSV; auto‑selects isoform FASTA from upstreams

## Create a New Run (TSV-based)

Once all stage configurations have been created : 

1) Copy and edit `config/templates/inputs.toml`:
   - Set `test_set_id`.
   - Fill `[run]` fields: `version`, `conda_env`, `data_dir`, `reads_file`, `genome_fa` (if needed), and shared inputs like `gtf`, `regions_tsv`, `junctions`, and TSS/TES BEDs.
   - Add stage-independent QC or defaults under `[qc]` if desired.

2) Edit stage templates as needed:
   - Add flags under `[run.stages.flags]`; leave blank for FLAIR defaults.
   - Do not set output `-o` flags; outputs are handled by the harness.

3) Build a TSV configuration file listing one test per line: first the base inputs, then the stage templates in execution order. Paths are relative to the TSV.

```text
config/templates/inputs.toml	config/templates/align.toml	config/templates/correct.toml	config/templates/collapse.toml
```

Run it with:

```bash
flair-test-suite path/to/cases.tsv
```

Outputs write to `./outputs/<test_set_id>/` and a `run_summary.log` is created there.

### Common Pitfalls & Fixes

- Missing `test_set_id`: add it at the top of `inputs.toml` (or under `[run]`).
- Wrong path base: remember stage TOML paths resolve relative to the TSV; data files resolve relative to `run.data_dir`.
- Stage order issues: ensure the TSV lists stages in valid order (e.g., `align` → `correct` → `collapse`).
- Conda env mismatch: set `run.conda_env` to an environment that has the appropriate FLAIR version and optional SQANTI tools.
- Permission problems: ensure you can write to `./outputs/` where you launch the command.

## Stage Ordering Rules

Order in the TSV is validated. Constraints:
- `regionalize` must follow `align`.
- `correct` must follow `align` (optionally after `regionalize`) and cannot follow `collapse`, `transcriptome`, `combine`, or `quantify`.
- `collapse` must follow `correct`.
- `transcriptome` must follow `align` and cannot follow `correct`.
- `combine` must follow `collapse` or `transcriptome` unless `flags.manifest` is provided.
- `quantify` must follow `combine`, `transcriptome`, or `collapse`.

## Output Layout & Caching

- Outputs live under `./outputs/<test_set_id>/<stage>/<signature>/`.
- Signatures incorporate tool version, flags, and hashes of input files/paths. Changing any of these yields a new signature directory.
- Re‑run logic:
  - If marker, primary output, and QC exist → stage is skipped.
  - If primary exists but QC is missing → QC is regenerated.
  - Otherwise → stage runs normally.

## Config Cookbook

- Targeted region runs: add `regionalize.toml` after `align.toml`; set `run.regions_tsv` (3 columns: `chr  start  end`).
- Multiple reads: use list or a comma‑string in `run.reads_file`.
- TED for collapse/transcriptome runs automatically; no config block is needed. Provide peak BEDs under `[run]` to enable precision/recall metrics; otherwise those metrics are `None`.

- SQANTI runs automatically when a conda env is available. Prefer setting `[run].sqanti_env` and optional `sqanti_cpus`.

- Combine with optional user manifest: set `flags.manifest = "my_manifest.tsv"`; upstream isoforms from this run are auto‑appended and deduplicated.
- Quantify with reads manifest: set `flags.manifest = "reads_manifest.tsv"`; isoform FASTA is auto‑selected from upstream outputs.

### Build the TSV (examples)

Place one workflow per line in your TSV, with base inputs first, then stage templates in execution order. See `config/templates/test_case_cli_template.tsv` for examples.

### Optional test-case notes

You can include commented context in `inputs.toml` explaining the test rationale (organism, samples, assembly, platform, region selection, etc.). This aids reproducibility and review.

## Tips

- Stage ordering is the TSV order.
- All outputs are forced to `./outputs/` by the harness regardless of `work_dir` in configs.

## Regionalized Effects

When you include a `regionalize` stage after `align`, downstream stages automatically switch to regionalized mode and operate per region. Key impacts:

- Inputs discovered: `correct` and `collapse`/`transcriptome` discover per‑region upstream files using the region index (`qc/regionalize/region_details.tsv`). Missing/empty region inputs are skipped with a warning.
- Filenames and tags: each region is identified by `{chrom}_{start}_{end}`. Downstream outputs use this tag:
  - `correct`: `{tag}_all_corrected.bed` plus `{tag}_all_inconsistent.bed`.
  - `collapse`/`transcriptome`: `{tag}.isoforms.bed` and `{tag}.isoforms.gtf` per region.
- QC outputs:
  - `correct` writes per‑region QC under `qc/<tag>/correct_qc.tsv` and `qc/<tag>/correct_splice_site_motifs.json`, plus an aggregate `qc/correct_qc.tsv` at the stage root with totals.
  - TED runs automatically and writes one row per region to `qc/ted/TED.tsv` and a transcriptome browser mapping at `qc/ted/transcriptome_browser/region_map.json` (used to open plots by region).
- Combine behavior: when `combine` runs, it auto‑discovers all `{tag}.isoforms.bed` files from upstream `collapse`/`transcriptome` and appends them to the manifest (alongside any user‑provided entries, deduplicated by path).
- Quantify behavior: `quantify` selects the correct isoform FASTA from upstream (`combine` preferred, otherwise `transcriptome`/`collapse`) with no change required for regionalized runs.
- Caching & signatures: including `regionalize` changes downstream signatures; a different set of inputs and the regionalize signature are hashed, so regionalized vs non‑regionalized runs won’t collide.
- Performance notes: each downstream stage runs once per region; expect more, smaller commands rather than one large run. Logs and summaries list region tags and counts.

Tip: keep region lists small for quick iteration; you can expand `regions_tsv` as your test case matures.
