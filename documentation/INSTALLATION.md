# TAUSO Installation Guide

The complete guide to installing TAUSO for all use cases.

---

## Quick Start (Most Users)

```bash
# One command - installs everything you need
mamba env create -f environment-dev.yml
mamba activate tauso

# Verify it worked
python -c "import tauso; print('Installed:', tauso.__file__)"
```

**That's it!** This installs TAUSO with all dependencies (runtime + development tools). Continue to [Data Setup](DATA_SETUP.md).

---

## Installation Methods

Choose based on your use case:

| Use Case | Method | Command |
|----------|--------|---------|
| **Development** (most users) | environment-dev.yml | `mamba env create -f environment-dev.yml` |
| **Production** (minimal runtime) | environment.yml + pip | `mamba env create -f environment.yml` → `pip install -e .` |
| **Reproducible research** | Lock file | `mamba create -n tauso --file conda-linux-64-dev.lock` → `pip install -e ".[dev]"` |

---

## Method 1: Development Setup (Recommended)

**Best for:** Daily use, contributing, running notebooks, training models

```bash
# Create environment with all dependencies
mamba env create -f environment-dev.yml

# Activate
mamba activate tauso

# Done! TAUSO is already installed with dev tools
```

**What you get:**
- ✅ TAUSO installed in editable mode
- ✅ All runtime dependencies (numpy, pandas, biopython, etc.)
- ✅ Development tools (pytest, ruff, jupyterlab)
- ✅ ML/analysis tools (matplotlib, seaborn, optuna, shap)
- ✅ Fast installer (uv)

**Next steps:** [Data Setup](DATA_SETUP.md)

---

## Method 2: Production/Minimal Runtime

**Best for:** Deployment, Docker containers, minimal installs

```bash
# Create minimal conda environment
mamba env create -f environment.yml
mamba activate tauso

# Install TAUSO (runtime dependencies only)
pip install -e .
```

**What you get:**
- ✅ TAUSO installed in editable mode
- ✅ Only runtime dependencies (smaller footprint)
- ❌ No pytest, matplotlib, or other dev tools

**To add dev tools later:**
```bash
pip install -e ".[dev]"
```

---

## Method 3: Reproducible Research

**Best for:** Published papers, exact version matching, scientific reproducibility

```bash
# Create from lock file (exact versions)
mamba create -n tauso --file conda-linux-64-dev.lock
mamba activate tauso

# Install TAUSO
pip install -e ".[dev]"
```

**What you get:**
- ✅ Exact package versions (bit-for-bit reproducible)
- ✅ Verified hashes
- ✅ Can recreate environment years later
- ⚠️ Linux-only (macOS needs separate lock file)

**For papers:** Document which lock file and git commit:
```
Environment: conda-linux-64-dev.lock from commit abc123
Repository: https://github.com/RedPenguin100/TAUSO
```

---

## Understanding `[dev]` Extras

TAUSO has optional development dependencies in `pyproject.toml`:

```bash
pip install -e .         # Runtime only (numpy, pandas, click, etc.)
pip install -e ".[dev]"  # Runtime + dev (pytest, matplotlib, optuna, etc.)
```

**[dev] extras include:**
- pytest, pytest-regressions, pytest-xdist (testing)
- matplotlib, seaborn (plotting)
- optuna (hyperparameter optimization)
- shap (model interpretability)
- pandas-stubs, scipy-stubs (type checking)

**Note:** `environment-dev.yml` automatically includes `[dev]` - you don't need to install separately!

---

## Setting Custom Data Directory

By default, TAUSO stores data in `~/.local/share/tauso`. To use a custom location:

### Option A: Environment-Specific (Recommended)

Sets variable only when conda environment is active:

```bash
# Activate environment
mamba activate tauso

# Set custom data directory for this environment
mamba env config vars set TAUSO_DATA_DIR=/path/to/your/data

# Reactivate to apply
mamba deactivate
mamba activate tauso

# Verify
echo $TAUSO_DATA_DIR
python -c "from tauso.data.data import get_data_dir; print(get_data_dir())"
```

**Benefits:**
- ✅ Only active in this environment
- ✅ Different data dirs for different envs
- ✅ Persists across sessions
- ✅ Clean and isolated

**Manage variables:**
```bash
mamba env config vars list                    # View all
mamba env config vars unset TAUSO_DATA_DIR    # Remove
```

### Option B: Shell-Wide

Sets variable globally in your shell:

```bash
# For bash
echo 'export TAUSO_DATA_DIR=/path/to/your/data' >> ~/.bashrc
source ~/.bashrc

# For zsh
echo 'export TAUSO_DATA_DIR=/path/to/your/data' >> ~/.zshrc
source ~/.zshrc
```

---

## Using `uv` (Faster Installer)

`uv` is a fast Python package installer (10-100x faster than pip). It's included in `environment-dev.yml`.

```bash
# Standard pip
pip install -e ".[dev]"

# Faster uv (drop-in replacement)
uv pip install -e ".[dev]"
```

Both produce identical results. Use whichever you prefer.

---

## Prerequisites

### Install Mamba (First-Time Setup)

If you don't have conda/mamba:

```bash
# Download Miniforge (includes mamba)
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"

# Install
bash Miniforge3-Linux-x86_64.sh

# Follow prompts, then restart shell
```

### System Dependencies

For `raccess` (RNA folding tool):

```bash
# Ubuntu/Debian
sudo apt-get install zlib1g-dev build-essential

# RHEL/CentOS
sudo yum install zlib-devel gcc gcc-c++
```

---

## Verification

After installation, verify everything works:

```bash
# Check TAUSO is installed
python -c "import tauso; print('TAUSO:', tauso.__file__)"

# Check CLI works
tauso --help

# Check data directory
python -c "from tauso.data.data import get_data_dir; print('Data dir:', get_data_dir())"

# If you installed [dev], check pytest
pytest --version
```

---

## Troubleshooting

### Import Error: No module named 'tauso'

You need to install TAUSO:
```bash
pip install -e ".[dev]"
```

### ModuleNotFoundError: No module named 'pytest'

You installed without `[dev]` extras:
```bash
pip install -e ".[dev]"
```

### environment-dev.yml installs but import still fails

Make sure you activated the environment:
```bash
mamba activate tauso
```

### Different data directory not working

Reactivate after setting the variable:
```bash
mamba env config vars set TAUSO_DATA_DIR=/path
mamba deactivate && mamba activate tauso
```

---

## Comparison Table

| Feature | environment-dev.yml | environment.yml + pip | Lock file |
|---------|-------------------|---------------------|-----------|
| **Commands** | 1 command | 2 commands | 2 commands |
| **TAUSO included** | ✅ Auto | ❌ Manual | ❌ Manual |
| **Dev tools** | ✅ Yes | ⚠️ With [dev] | ⚠️ With [dev] |
| **Reproducibility** | ⚠️ Approximate | ⚠️ Approximate | ✅ Exact |
| **Cross-platform** | ✅ Yes | ✅ Yes | ❌ Linux-only |
| **Install time** | Fast | Fast | Fastest |
| **Disk usage** | ~2 GB | ~1 GB | ~2 GB |
| **Best for** | Development | Production | Research papers |

---

## What's Next?

After installation:

1. **Set up data** → [Data Setup Guide](DATA_SETUP.md)
2. **Quick reference** → [Setup Reference](SETUP_REFERENCE.md)
3. **Start using TAUSO** → See main [README](../README.md)

---

## Advanced: Production Deployment

For HPC clusters or reproducible production pipelines, see [Reproducible Setup Guide](feature_run/ReproducibleSetup.md).
