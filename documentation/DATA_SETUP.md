# TAUSO Data Setup Guide

Download and prepare datasets for TAUSO feature extraction and model training.

**Prerequisites:** TAUSO must be installed → [Installation Guide](INSTALLATION.md)

---

## Custom Data Directory (Optional)

**Default:** `~/.local/share/tauso`

To use a custom location:
```bash
mamba env config vars set TAUSO_DATA_DIR=/path/to/your/data
mamba deactivate && mamba activate tauso
```

**→ Details:** [Setting Data Directory](INSTALLATION.md#setting-custom-data-directory)

---

## Installation Approaches

Choose the approach that fits your needs:

### Option A: One Command (Fastest)

Best for: First-time users, quick setup

```bash
# Set up everything in one go
tauso setup-all

# Build cell-line context files
tauso build-cell-context

# (Optional) Training data
python -m notebooks.data.OligoAI.assign_canonical_gene
python -m notebooks.utils.data
```

**Time:** 1-2 hours (mostly bowtie indexing)

---

### Option B: Incremental (Recommended for Custom Setups)

Best for: Custom data directory, step-by-step verification, selective features

```bash
# 1. Set custom data directory (optional)
export TAUSO_DATA_DIR=/path/to/your/data

# 2. Core genome data (required)
tauso setup-genome          # Downloads genome + annotations (5-10 min)

# 3. Off-target support (optional, but recommended)
tauso setup-bowtie          # Builds alignment index (30-60 min)

# 4. RNA folding (required for structure features)
tauso setup-raccess         # Compiles folding tool (2-5 min)

# 5. Omics data (required for cell-context features)
tauso setup-depmap          # Cell line expression (5-10 min)
tauso setup-mrna-halflife   # mRNA stability (1-2 min)
tauso setup-tgcn            # tRNA gene counts (<1 min)
tauso setup-attract         # RBP motifs (1-2 min)
tauso setup-riboseq         # Ribosome profiling (2-5 min)

# 6. Build cell context (required for per-cell features)
tauso build-cell-context    # Expression + CAI weights (5-10 min)

# 7. Training data (only if training models)
python -m notebooks.data.OligoAI.assign_canonical_gene
python -m notebooks.utils.data
```

**Benefits of incremental approach:**
- Verify each step completes successfully
- Skip optional components you don't need
- Use custom data directory more easily
- Resume if interrupted

---

## Step-by-Step Breakdown

If you want more control or need to troubleshoot, here's what each command does:

### 1. Core Reference Data (Required)

#### Genome & Annotations
```bash
tauso setup-genome [--genome GRCh38]
```
- Downloads FASTA genome sequence and GTF/GFF annotations
- Builds gffutils database for fast gene queries
- **Time:** 5-10 minutes
- **Disk:** ~3 GB (GRCh38)
- **Supported genomes:** GRCh38 (human), GRCm39 (mouse)

#### Bowtie Index
```bash
tauso setup-bowtie [-t THREADS] [--mem-per-thread MB]
```
- Creates bowtie alignment index for off-target search
- **Time:** 30-60 minutes (single-threaded), scales with threads
- **Disk:** ~4 GB
- **Tip:** Use `-t $(nproc)` for faster indexing on multi-core systems
- **Memory:** Default 800 MB/thread, increase if you have RAM

#### RNA Folding Tool
```bash
tauso setup-raccess
```
- Installs raccess (RNA accessibility prediction)
- **Required:** `zlib1g-dev` or equivalent must be installed first
- **Time:** 2-5 minutes
- **Note:** Distribution restricted by license, must compile locally

---

### 2. Omics Data (Required for Feature Calculation)

You can install all omics datasets at once or individually:

#### All Omics (Recommended)
```bash
tauso setup-omics [--force]
```
Includes: DepMap, mRNA half-life, tGCN (tAI), ATtRACT RBP motifs, ribo-seq

#### Individual Omics Datasets

**DepMap (Cell Line Expression)**
```bash
tauso setup-depmap
```
- Downloads DepMap Public 25Q3 from Zenodo
- Cell line metadata, profiles, and expression data
- Converts CSV to Parquet for fast loading
- **Disk:** ~500 MB (Parquet)

**mRNA Half-Life**
```bash
tauso setup-mrna-halflife
```
- TTDB (Transcriptome Turnover Database) data
- **Disk:** ~50 MB

**tRNA Gene Copy Numbers (tGCN)**
```bash
tauso setup-tgcn
```
- Fetched from GtRNAdb
- Used for tAI (tRNA Adaptation Index) features
- **Disk:** <1 MB

**ATtRACT RBP Motifs**
```bash
tauso setup-attract
```
- RNA-binding protein motif database
- **Disk:** ~10 MB

**Ribo-Seq BigWig**
```bash
tauso setup-riboseq
```
- 40S ribosome profiling data
- **Disk:** ~200 MB

---

### 3. Cell Context (Required for Per-Cell Features)

```bash
tauso build-cell-context [CELL_NAMES...] [--reset] [--no-cai]
```

This command:
1. Registers cell lines in a cohort (uses 30 default cell lines if none specified)
2. Extracts per-cell expression files from DepMap
3. Computes CAI (Codon Adaptation Index) weights per cell line
4. Sets up tGCN data

**Examples:**
```bash
# Use default cohort (30 cell lines)
tauso build-cell-context

# Custom cell lines (append to existing)
tauso build-cell-context HEPG2 HELA A549

# Replace cohort
tauso build-cell-context --reset HEPG2 HELA A549

# Skip CAI computation
tauso build-cell-context --no-cai
```

**Output files:**
- `cell_cohort.json` - List of cell lines in your cohort
- `processed_expression/*.csv` - Per-cell expression files
- `cai_weights.json` - CAI weights for each cell line

---

### 4. Training Data (Only if Training Models)

If you plan to train ML models or run tests, prepare the OligoAI training data:

```bash
# Assign canonical genes to oligonucleotides
python -m notebooks.data.OligoAI.assign_canonical_gene

# Create index for fast lookups
python -m notebooks.utils.data
```

**Time:** 1-2 minutes

---

## Command Comparison

| Command | What It Does | Required For | Time |
|---------|--------------|-------------|------|
| `setup-all` | Genome + bowtie + omics + raccess | Everything | 60-90 min |
| `setup-genome` | Download genome FASTA + annotations | Core features | 5-10 min |
| `setup-bowtie` | Index genome for alignment | Off-target analysis | 30-60 min |
| `setup-raccess` | Install RNA folding tool | Folding features | 2-5 min |
| `setup-omics` | All omics datasets | Cell-context features | 10-20 min |
| `setup-depmap` | Cell line expression | Expression features | 5-10 min |
| `setup-mrna-halflife` | mRNA stability data | Half-life features | 1-2 min |
| `setup-tgcn` | tRNA gene copy numbers | tAI features | <1 min |
| `setup-attract` | RBP motif database | RBP features | 1-2 min |
| `setup-riboseq` | Ribosome profiling | Translation features | 2-5 min |
| `build-cell-context` | Cell cohort + expression + CAI | Per-cell features | 5-10 min |

---

## Verification

Check that everything is set up correctly:

```bash
# Check current data directory
python -c "from tauso.data.data import get_data_dir; print(get_data_dir())"
# Default: ~/.local/share/tauso
# Or: your custom TAUSO_DATA_DIR if set

# List downloaded files
ls -lh $(python -c "from tauso.data.data import get_data_dir; print(get_data_dir())")

# Expected files after full setup:
# - GRCh38.fa (genome FASTA)
# - GRCh38.gtf.db (annotation database)
# - GRCh38_bowtie_index/ (bowtie index files)
# - raccess/ (RNA folding tool)
# - Model.csv, OmicsProfiles.csv (DepMap)
# - OmicsExpressionTPMLogp1HumanAllGenesStranded.parquet
# - mrna_half_life.parquet
# - human_tgcn_hsapi38.csv
# - attract/ (RBP motifs)
# - human_unselected_40S.RiboProElong.bw
# - processed_expression/ (per-cell files)
# - cell_cohort.json
# - cai_weights.json
```

---

## Troubleshooting

### Bowtie indexing fails with memory error
Reduce memory per thread:
```bash
tauso setup-bowtie --mem-per-thread 400
```

### raccess compilation fails
Install development libraries:
```bash
# Ubuntu/Debian
sudo apt-get install zlib1g-dev build-essential

# RHEL/CentOS
sudo yum install zlib-devel gcc gcc-c++
```

### DepMap download is slow
The Zenodo mirror should be fast. If not, check your internet connection or try again later.

### Missing cell line in cohort
Check available cell lines:
```bash
# Cell lines must exist in OmicsProfiles.csv
# Search is fuzzy and case-insensitive
tauso add-cell "HepG2" "HeLa"  # Will find HEPG2, HELA
```

### Data directory location
**Default:** `~/.local/share/tauso`

**Change temporarily:**
```bash
export TAUSO_DATA_DIR=/path/to/data
tauso setup-all  # Will use the new location
```

**Change permanently:** See "Setting Data Directory" section at the top of this guide.

---

## What's Next?

After data setup is complete, you can:

1. **Run feature extraction** on your ASO designs
2. **Analyze off-targets:** `tauso run-off-target SEQUENCE`
3. **Train ML models** (see `notebooks/models/`)
4. **Design ASOs** using the `design_asos` function

See the main README for usage examples.
