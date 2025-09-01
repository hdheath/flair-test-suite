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

- Have FLAIR installed in a conda env you will use (e.g., `flair`).
- Place all input data for a run under a single `data_dir` (you can symlink large files to avoid duplication).

  ```bash
  ln -s /path/to/shared/genome.fa /path/to/data_dir/genome.fa
  ```

---

## Running the Test Suite

Define runs with a single TSV named `flair_test_suite_config.tsv` (two sections) per test set. The first section lists inputs as `key<TAB>value`; the second lists stages `stage<TAB>flags`. A new case starts at each `align` row. Paths resolve relative to the TSV file.

1) Copy and edit the template into your config file:
   - `cp config/templates/combined_cases.tsv flair_test_suite_config.tsv`
   - In the `key\tvalue` section, set `test_set_id`, `version`, `flair_env` (Conda env containing FLAIR), `data_dir`, `reads_file`, and shared inputs (`genome_fa`, `gtf`, optional `junctions`, `sqanti_env`).
   - In the `stage\tflags` section, list stages in order. Flags are written exactly as in the FLAIR CLI, but comma-separated within the row. Both `--opt value` and `--opt=value` forms are accepted (e.g., `--nvrna, --threads 8`). Use `region_test` and pass `--regions-tsv regions.tsv` for region runs.

2) Run the suite:

```bash
flair-test-suite flair_test_suite_config.tsv
```

Notes:
- Outputs write to `./outputs/<test_set_id>/` with a `run_summary.log` per test set.
- Region-scoped runs: add a `region_test` row after `align` and pass `--regions-tsv=<path>`; region TSV is 3 columns: `chr  start  end`.

For more configuration details and examples, see: `docs/configurations.md`.

---

## Contributing & License

* Bug reports and PRs welcome — see **`CONTRIBUTING.md`**.
* Code of Conduct: **`CODE_OF_CONDUCT.md`**.
* © 2025 **Harrison Heath / Brooks Lab** – released under the **MIT License** (`LICENSE`).

Further details: [FLAIR Test Suite Overview](docs/overview.md)

Happy testing 🚀
