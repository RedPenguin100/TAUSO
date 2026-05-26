# TAUSO Setup Quick Reference

## Setting Custom Data Directory

### Method 1: Environment-Specific (Recommended)
```bash
# Activate environment
mamba activate tauso

# Set variable for this environment only
mamba env config vars set TAUSO_DATA_DIR=/path/to/your/data

# Reactivate to apply
mamba deactivate && mamba activate tauso

# Verify
echo $TAUSO_DATA_DIR
```

### Method 2: Shell-Wide
```bash
# For bash
echo 'export TAUSO_DATA_DIR=/path/to/your/data' >> ~/.bashrc
source ~/.bashrc

# For zsh
echo 'export TAUSO_DATA_DIR=/path/to/your/data' >> ~/.zshrc
source ~/.zshrc
```

### Method 3: Temporary (Current Session)
```bash
export TAUSO_DATA_DIR=/path/to/your/data
```

### Check Current Location
```bash
python -c "from tauso.data.data import get_data_dir; print(get_data_dir())"
```

**Default location:** `~/.local/share/tauso`

### Manage Environment Variables
```bash
# List all env vars for current environment
mamba env config vars list

# Remove variable
mamba env config vars unset TAUSO_DATA_DIR
```

---

## Setup Commands

### All-in-One (Recommended for First-Time Users)
```bash
tauso setup-all                    # Everything: genome + bowtie + omics + raccess
tauso build-cell-context           # Cell context: cohort + expression + CAI
```

### Incremental (Recommended for Custom Setups)

#### 1. Core Genome Data (Required)
```bash
tauso setup-genome [--genome GRCh38]     # 5-10 min, ~3 GB
tauso setup-bowtie [-t THREADS]          # 30-60 min, ~4 GB
tauso setup-raccess                      # 2-5 min, ~50 MB
```

#### 2. Omics Data (Required for Cell Features)
```bash
# All omics at once
tauso setup-omics                        # 10-20 min, ~1 GB

# Or individually:
tauso setup-depmap                       # Cell line expression (5-10 min, ~500 MB)
tauso setup-mrna-halflife                # mRNA stability (1-2 min, ~50 MB)
tauso setup-tgcn                         # tRNA gene counts (<1 min, <1 MB)
tauso setup-attract                      # RBP motifs (1-2 min, ~10 MB)
tauso setup-riboseq                      # Ribosome profiling (2-5 min, ~200 MB)
```

#### 3. Cell Context (Required for Per-Cell Features)
```bash
tauso build-cell-context                 # 5-10 min
# Or with custom cells:
tauso build-cell-context HEPG2 HELA A549
tauso build-cell-context --reset HEPG2 HELA  # Replace existing cohort
```

#### 4. Training Data (Only if Training Models)
```bash
python -m notebooks.data.OligoAI.assign_canonical_gene
python -m notebooks.utils.data
```

---

## Command Options

### Common Flags
```bash
--force              # Force redownload/rebuild
--help               # Show help for any command
```

### setup-genome
```bash
--genome GRCh38      # Genome version (GRCh38, GRCm39)
--force              # Force redownload
--remove-gz          # Remove compressed files after extraction
```

### setup-bowtie
```bash
-t, --threads N      # Number of CPU threads (default: 1)
--mem-per-thread MB  # Memory per thread (default: 800)
--force              # Force rebuild
```

### setup-all
```bash
--genome GRCh38      # Genome version
-t, --threads N      # Threads for bowtie
--mem-per-thread MB  # Memory per thread for bowtie
--force              # Force redownload all
```

### build-cell-context
```bash
CELL_NAMES...        # Cell line names to add (e.g., HEPG2 HELA)
--reset              # Clear existing cohort before adding
--no-cai             # Skip CAI weight computation
--genome GRCh38      # Genome version
--force              # Force CAI recomputation
```

---

## Verification Commands

```bash
# Check TAUSO version
python -c "import tauso; print(tauso.__file__)"

# Check data directory
python -c "from tauso.data.data import get_data_dir; print(get_data_dir())"

# List downloaded files
ls -lh $(python -c "from tauso.data.data import get_data_dir; print(get_data_dir())")

# Test genome setup
python -c "from tauso.data.data import load_gtf_db; db = load_gtf_db('GRCh38'); print('Genes:', db.count_features_of_type('gene'))"

# Test off-target search
tauso run-off-target ACGTACGTACGTACGTACGT
```

---

## Typical Workflows

### Minimal Setup (Just Genome + Folding)
```bash
export TAUSO_DATA_DIR=/data/tauso
tauso setup-genome
tauso setup-raccess
# Can now run basic features (no off-target, no cell-specific)
```

### Off-Target Analysis Setup
```bash
export TAUSO_DATA_DIR=/data/tauso
tauso setup-genome
tauso setup-bowtie -t $(nproc)  # Use all CPU cores
# Can now run: tauso run-off-target SEQUENCE
```

### Full Feature Extraction Setup
```bash
export TAUSO_DATA_DIR=/data/tauso
tauso setup-all -t $(nproc)
tauso build-cell-context
# Can now calculate all features
```

### ML Training Setup
```bash
export TAUSO_DATA_DIR=/data/tauso
tauso setup-all -t $(nproc)
tauso build-cell-context
python -m notebooks.data.OligoAI.assign_canonical_gene
python -m notebooks.utils.data
# Can now train models
```

---

## Disk Space Requirements

| Component | Size | Required For |
|-----------|------|--------------|
| Genome (GRCh38) | ~3 GB | Everything |
| Bowtie index | ~4 GB | Off-target analysis |
| DepMap | ~500 MB | Cell expression features |
| mRNA half-life | ~50 MB | Stability features |
| ATtRACT RBP | ~10 MB | RBP features |
| Ribo-seq | ~200 MB | Translation features |
| tGCN | <1 MB | tAI features |
| raccess | ~50 MB | Folding features |
| Processed expression | ~100 MB | Per-cell features |
| CAI weights | ~5 MB | CAI features |
| **Total (full)** | **~8-10 GB** | All features |

---

## Time Estimates (Typical Hardware)

| Command | Time | Bottleneck |
|---------|------|------------|
| setup-genome | 5-10 min | Download speed |
| setup-bowtie (1 thread) | 30-60 min | CPU |
| setup-bowtie (16 threads) | 5-10 min | CPU |
| setup-raccess | 2-5 min | Compilation |
| setup-depmap | 5-10 min | Download + conversion |
| setup-mrna-halflife | 1-2 min | Download |
| setup-tgcn | <1 min | Network fetch |
| setup-attract | 1-2 min | Download |
| setup-riboseq | 2-5 min | Download |
| build-cell-context | 5-10 min | CAI computation |
| **setup-all (16 threads)** | **30-45 min** | Bowtie indexing |
| **Full setup** | **45-60 min** | - |

---

## Troubleshooting Quick Fixes

### raccess won't compile
```bash
sudo apt-get install zlib1g-dev build-essential
```

### Bowtie out of memory
```bash
tauso setup-bowtie --mem-per-thread 400  # Reduce from 800
```

### Cell line not found
```bash
# Search is case-insensitive and fuzzy
tauso add-cell "HepG2"    # Finds HEPG2
tauso add-cell "HeLa"     # Finds HELA
```

### Wrong data directory
```bash
# Check where it's actually looking
python -c "from tauso.data.data import get_data_dir; print(get_data_dir())"

# Make sure TAUSO_DATA_DIR is set correctly
echo $TAUSO_DATA_DIR
```

### Fresh start
```bash
# Remove everything and start over
rm -rf $(python -c "from tauso.data.data import get_data_dir; print(get_data_dir())")
tauso setup-all --force
```

---

## See Also

- **Complete guide:** [DATA_SETUP.md](DATA_SETUP.md)
- **Quick start:** [QUICKSTART.md](QUICKSTART.md)
- **Installation:** [INSTALLATION.md](INSTALLATION.md)
- **Production setup:** [feature_run/ReproducibleSetup.md](feature_run/ReproducibleSetup.md)
