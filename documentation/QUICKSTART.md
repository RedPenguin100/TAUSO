# TAUSO Quick Start

The short path from a clean machine to a working TAUSO. For the reasoning behind each step, follow the linked guides.

## 1. Install

```bash
mamba env create -f environment-dev.yml
mamba activate tauso
python -c "import tauso; print(tauso.__file__)"
```

Building the environment downloads a large dependency set and takes several minutes. See [Installation](INSTALLATION.md) for the runtime-only and reproducible variants.

## 2. Choose a data directory

A full data setup is on the order of 10 GB. Decide where it lives **before** step 3, or it goes to `~/.local/share/tauso`:

```bash
mamba env config vars set TAUSO_DATA_DIR=/path/to/your/data
mamba deactivate && mamba activate tauso
```

## 3. Download and build the data

```bash
tauso setup-all            # genome + bowtie off-target index + omics + raccess
tauso build-cell-context   # per-cell expression, CAI weights, human tGCN
```

Most of the time here is the genome download and the off-target index build — expect the better part of an hour. Pass `-t $(nproc)` to `setup-all` to parallelize indexing. See [Data Setup](DATA_SETUP.md) for the step-by-step path and what each dataset is for.

## Next

- [Data Setup](DATA_SETUP.md) — install datasets selectively and verify them.
- [Training Data](TRAINING_DATA.md) — prepare the dataset for model training.
- [Setup Reference](SETUP_REFERENCE.md) — command, disk, and time cheat sheet.
