# Installing TAUSO

TAUSO is distributed as a conda/pip package. Which environment you build depends on what you intend to do:

- **Develop and run the models** — the full environment with development and analysis tooling. What most contributors want.
- **Use and generate ASOs** — a minimal runtime environment.
- **Reproduce a published result** — a pinned lock file.

Installation supports both `pip` and `uv`; the commands below use `pip`, and `uv pip ...` is a drop-in substitute.

## Prerequisites

A conda-compatible package manager (we use `mamba`). If you don't have one, install Miniforge:

```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh
# restart your shell afterwards
```

`raccess`, the RNA-folding tool used during data setup, is compiled from source and needs a C/C++ toolchain and zlib headers:

```bash
# Debian/Ubuntu
sudo apt-get install build-essential zlib1g-dev
# RHEL/CentOS
sudo yum install gcc gcc-c++ zlib-devel
```

## Develop and run the models

`environment-dev.yml` builds the conda environment and installs TAUSO in editable mode with the `[dev]` extras in a single step. Resolving and downloading the full dependency set takes several minutes, and longer on a cold conda cache.

```bash
mamba env create -f environment-dev.yml
mamba activate tauso
```

Verify:

```bash
python -c "import tauso; print(tauso.__file__)"
tauso --help
pytest --version
```

## Use and generate ASOs

For the library at run time only, build the minimal environment and install TAUSO without the dev extras:

```bash
mamba env create -f environment.yml
mamba activate tauso
pip install -e .
```

The development extras can be added later with `pip install -e ".[dev]"`.

## Reproduce a published result

The lock file pins exact package versions (with hashes) so the environment can be rebuilt identically — for example, years later from a paper's reported commit. The lock files are Linux-only; macOS needs its own.

```bash
mamba create -n tauso --file conda-linux-64-dev.lock
mamba activate tauso
pip install -e ".[dev]"
```

When reporting an environment in a paper, record both the lock file and the git commit it came from.

## The `[dev]` extras

The `[dev]` optional dependencies add the test suite and the analysis/plotting stack on top of the runtime requirements — for example `pytest`, `matplotlib`, and `pandas-stubs`. The full, current list is in `pyproject.toml`. `environment-dev.yml` already includes them, so a development environment needs no separate step.

## Choosing a data directory

TAUSO stores downloaded genomes and datasets in `~/.local/share/tauso` by default. A full setup is on the order of 10 GB (see [Data Setup](DATA_SETUP.md)), so decide where this should live **before** the first `tauso setup-all` — otherwise it populates the default location. To put it elsewhere, set `TAUSO_DATA_DIR`.

Per environment (scoped to the `tauso` environment, persists across sessions):

```bash
mamba activate tauso
mamba env config vars set TAUSO_DATA_DIR=/path/to/your/data
mamba deactivate && mamba activate tauso        # reactivate to apply
```

Shell-wide:

```bash
echo 'export TAUSO_DATA_DIR=/path/to/your/data' >> ~/.bashrc   # or ~/.zshrc
```

Confirm the active location:

```bash
python -c "from tauso.data.data import get_data_dir; print(get_data_dir())"
```

## Next

- [Data Setup](DATA_SETUP.md) — download genomes and datasets.
- [Setup Reference](SETUP_REFERENCE.md) — command, disk, and time cheat sheet.
- [Reproducible production setup](feature_run/ReproducibleSetup.md) — HPC and pipeline deployment.
