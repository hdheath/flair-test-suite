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

* Users must already have downloaded the FLAIR conda env they would like to use.
* All input data used for a run must be located in the same directory.

  * Tip: you can symlink large files to the input directory to avoid duplication.

    ```bash
    # pseudo command to create a symbolic link instead of copying large files
    ln -s /path/to/shared/genome.fa /path/to/data_dir/genome.fa
    ```
* Each configuration file must be created by the user.

For help creating these, visit : 
[Guide to Creating Configurations for Test Cases](docs/configurations.md)

---

## Running the Test Suite

You can run individual test-cases or batch multiple configurations.

### Single config

```bash
flair-test-suite path/to/your_config.toml
```

This will execute the workflow defined in `your_config.toml` and place results under `outputs/<run_id>/`.

By default the run summary logs paths relative to the working directory. Use
`--absolute-paths` if you prefer full paths:

```bash
flair-test-suite --absolute-paths path/to/your_config.toml
```

### Batch runs from a list of configs

If you have many configs, create a text file listing each TOML path (one per line):

```text
configs/demo.toml
configs/sample1.toml
configs/sample2.toml
```

Then launch:

```bash
flair-test-suite path/to/runs.txt
```

---

## Contributing & License

* Bug reports and PRs welcome — see **`CONTRIBUTING.md`**.
* Code of Conduct: **`CODE_OF_CONDUCT.md`**.
* © 2025 **Harrison Heath / Brooks Lab** – released under the **MIT License** (`LICENSE`).

Further details: [FLAIR Test Suite Overview](docs/overview.md)

Happy testing 🚀

