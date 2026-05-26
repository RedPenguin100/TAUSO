# TAUSO Quick Start

Get TAUSO running in under 5 minutes (plus data download time).

---

## 1. Install (2 minutes)

```bash
# Create environment and install TAUSO
mamba env create -f environment-dev.yml
mamba activate tauso

# Verify
python -c "import tauso; print('✓ TAUSO installed')"
```

**Need help?** → [Installation Guide](INSTALLATION.md)

---

## 2. Set Data Directory (Optional)

```bash
# Set custom location (environment-specific)
mamba env config vars set TAUSO_DATA_DIR=/path/to/your/data
mamba deactivate && mamba activate tauso
```

**Default:** `~/.local/share/tauso`

---

## 3. Download Data (1-2 hours)

```bash
# One command - downloads everything
tauso setup-all

# Build cell context
tauso build-cell-context
```

**What this does:**
- Downloads human genome (GRCh38) + annotations (~3 GB)
- Builds bowtie index for off-target search (~4 GB)
- Downloads omics data (DepMap, mRNA half-life, etc.) (~1 GB)
- Installs raccess (RNA folding tool)
- Builds cell-line expression profiles

**Need step-by-step?** → [Data Setup Guide](DATA_SETUP.md)

---

## 4. Verify Setup

```bash
# Check data directory
python -c "from tauso.data.data import get_data_dir; print('Data:', get_data_dir())"

# Test off-target search
tauso run-off-target ACGTACGTACGTACGTACGT

# Should show genomic hits for this sequence
```

---

## 5. Start Using TAUSO

### Off-Target Analysis
```bash
tauso run-off-target ACGTACGTACGTACGTACGT -m 3 -o results.csv
```

### Python API
```python
from tauso import design_asos

# Design ASOs for a gene
asos = design_asos(
    gene_name="KRAS",
    genome="GRCh38",
    chemistry="MOE",
    include_details=True
)

print(asos.head())
```

---

## Troubleshooting

### "ModuleNotFoundError: No module named 'tauso'"
```bash
# Make sure environment is activated
mamba activate tauso
```

### "No such file: GRCh38.fa"
```bash
# Run data setup first
tauso setup-all
```

### Need more help?
- Installation issues → [INSTALLATION.md](INSTALLATION.md)
- Data setup issues → [DATA_SETUP.md](DATA_SETUP.md)
- Command reference → [SETUP_REFERENCE.md](SETUP_REFERENCE.md)

---

## What's Next?

- **Calculate features** → See `notebooks/features/` for examples
- **Train models** → See `notebooks/models/` for ML pipelines
- **API docs** → Check docstrings in `src/tauso/`
- **Production setup** → [Reproducible Setup](feature_run/ReproducibleSetup.md)
