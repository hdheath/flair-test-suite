# FLAIR Test Suite 🧪

*A config-driven regression and QC harness for*
[FLAIR](https://github.com/BrooksLabUCSC/flair) **Ver. ≥ 2.0**

---

## Table of Contents

1. [Installation & Setup](#installation--setup-single-step-conda--pip)
2. [Before Running the Test Suite](#before-running-the-test-suite)
3. [Running the Test Suite](#running-the-test-suite)
4. [Contributing & License](#contributing--license)

---

## Installation & Setup

> Requires Conda or Mamba.

```bash
git clone https://github.com/hdheath/flair-test-suite.git
cd flair-test-suite
conda env create -f flair-test-suite.yaml
conda activate flair-test-suite
```

### Verify

```bash
which flair-test-suite  # should print path inside env
```

---

## Before Running the Test Suite

* Users must already have downloaded FLAIR in a conda env they would like to use.
* All input data used for a run must be located in the same directory.

  * Tip: you can symlink large files to the input directory to avoid duplication.

    ```bash
    # pseudo command to create a symbolic link instead of copying large files
    ln -s /path/to/shared/genome.fa /path/to/data_dir/genome.fa
    ```

---

## Running the Test Suite

Runs are defined by a TSV configuration file that points to TOML templates. Paths are resolved relative to the TSV file.

1) Copy and edit `config/templates/inputs.toml` (base inputs).
   - Set `test_set_id`.
   - Fill `run.version`, `run.conda_env`, `run.data_dir`, `run.reads_file`, and any shared inputs (`gtf`, `regions_tsv`, `junctions`, TSS/TES BEDs).
2) Choose stage templates from `config/templates/` and edit flags as needed:
   - `align.toml`, `correct.toml`, optional `regionalize.toml`, then `collapse.toml` or `transcriptome.toml`.
   - Optional downstream: `combine.toml`, `quantify.toml` (see configuration templates in the same folder).

For help creating these, visit:
[Guide to Creating Configurations for Test Cases](docs/configurations.md)

3) Create a TSV configuration (one test case per line):

```text
config/templates/inputs.toml	config/templates/align.toml	config/templates/correct.toml	config/templates/collapse.toml
```

4) Run the suite:

```bash
flair-test-suite path/to/cases.tsv
```

Notes:
- Outputs always write to `./outputs/<test_set_id>/`. A `run_summary.log` is created per test set.
- For targeted region runs, add `regionalize.toml` and provide a 3-column region TSV (see `config/templates/region_template.tsv`).

---

## Contributing & License

* Bug reports and PRs welcome — see **`CONTRIBUTING.md`**.
* Code of Conduct: **`CODE_OF_CONDUCT.md`**.
* © 2025 **Harrison Heath / Brooks Lab** – released under the **MIT License** (`LICENSE`).

Further details: [FLAIR Test Suite Overview](docs/overview.md)

Happy testing 🚀
