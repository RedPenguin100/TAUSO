# TAUSO Setup Reference

A command, disk, and time cheat sheet. For explanations, see [Data Setup](DATA_SETUP.md).

## Data directory

```bash
# Per environment (persists, scoped to the tauso env)
mamba activate tauso
mamba env config vars set TAUSO_DATA_DIR=/path/to/your/data
mamba deactivate && mamba activate tauso

# Shell-wide
echo 'export TAUSO_DATA_DIR=/path/to/your/data' >> ~/.bashrc   # or ~/.zshrc

# Current session only
export TAUSO_DATA_DIR=/path/to/your/data

# Where is it now?
python -c "from tauso.data.data import get_data_dir; print(get_data_dir())"

# Manage env vars
mamba env config vars list
mamba env config vars unset TAUSO_DATA_DIR
```

Default location: `~/.local/share/tauso`. Set this before the first `tauso setup-all`.

## Setup commands

```bash
# Everything the package downloads/builds (see TODO in DATA_SETUP.md — omits the two below)
tauso setup-all                          # genome + bowtie + omics + raccess
tauso build-cell-context                 # cohort + per-cell expression + CAI + tGCN

# Genome (required)
tauso setup-genome [--genome GRCh38]     # genome FASTA + annotations + gffutils DB
tauso setup-bowtie [-t THREADS]          # off-target alignment index
tauso setup-raccess                      # RNA folding tool (compiled locally)

# Omics
tauso setup-omics                        # all of the below
tauso setup-depmap                       # cell-line expression
tauso setup-mrna-halflife                # mRNA stability
tauso setup-tgcn                         # tRNA gene copy numbers (for tAI)
tauso setup-attract                      # RBP motifs
tauso setup-riboseq                      # 40S scanning + 80S elongation (Wagner 2020)

# Cell context
tauso build-cell-context                 # default cohort
tauso build-cell-context HEPG2 HELA A549 # custom cohort (append)
tauso build-cell-context --reset HEPG2   # replace cohort

# Training data (only to train models / run tests — see TRAINING_DATA.md)
python -m notebooks.data.OligoAI.assign_canonical_gene
python -m notebooks.utils.data
python -m notebooks.data.OligoAI.process_data
```

## Command options

```bash
--force              # force redownload/rebuild
--help               # help for any command

# setup-genome
--genome GRCh38      # GRCh38 (human); GRCm39 (mouse) is WIP, not fully supported
--remove-gz          # delete compressed files after extraction

# setup-bowtie
-t, --threads N      # CPU threads (default 1)
--mem-per-thread MB  # memory per thread (default 800)

# setup-all
--genome GRCh38
-t, --threads N
--mem-per-thread MB

# build-cell-context
CELL_NAMES...        # cell lines to add (e.g. HEPG2 HELA)
--reset              # clear the cohort first
--no-cai             # skip CAI weights
```

## Verification

```bash
# Package
python -c "import tauso; print(tauso.__file__)"

# Data directory and contents
python -c "from tauso.data.data import get_data_dir; print(get_data_dir())"
ls -lh "$(python -c 'from tauso.data.data import get_data_dir; print(get_data_dir())')"

# Genome / annotation DB loads
python -c "from tauso.data.data import load_gtf_db; db = load_gtf_db('GRCh38'); print('genes:', db.count_features_of_type('gene'))"
```

## Disk space (GRCh38)

| Component | Size | Required for |
|-----------|------|--------------|
| Genome FASTA | ~3 GB | Everything |
| Annotation DB (gffutils) | ~3 GB | Everything |
| Bowtie index | ~3 GB | Off-target analysis |
| DepMap | ~1 GB | Cell-expression features |
| Ribo-seq (40S + 80S) | ~20 MB | Translation features |
| mRNA half-life | ~150 MB | Stability features |
| Processed expression | ~100 MB | Per-cell features |
| ATtRACT / raccess / tGCN / CAI | <100 MB | Respective features |
| **Total (full setup)** | **~10–12 GB** | All features |

The genome and the off-target index account for most of this — they are the cost of doing genome-wide off-target analysis.

## Time estimates (typical hardware)

| Step | Time | Bottleneck |
|------|------|------------|
| setup-genome | 5–10 min | Download |
| setup-bowtie (1 thread) | 30–60 min | CPU |
| setup-bowtie (16 threads) | 5–10 min | CPU |
| setup-raccess | 2–5 min | Compilation |
| setup-omics | 10–20 min | Download + conversion |
| build-cell-context | 5–10 min | CAI computation |
| **setup-all (16 threads)** | **~30–45 min** | Bowtie indexing |

## Typical workflows

```bash
# Minimal: genome + folding only (basic features; no off-target, no cell-specific)
tauso setup-genome
tauso setup-raccess

# Off-target analysis
tauso setup-genome
tauso setup-bowtie -t $(nproc)

# Full feature extraction
tauso setup-all -t $(nproc)
tauso build-cell-context

# Model training (adds the training data — see TRAINING_DATA.md)
tauso setup-all -t $(nproc)
tauso build-cell-context
python -m notebooks.data.OligoAI.assign_canonical_gene
python -m notebooks.utils.data
python -m notebooks.data.OligoAI.process_data
```

## See also

- [Installation](INSTALLATION.md)
- [Quick Start](QUICKSTART.md)
- [Data Setup](DATA_SETUP.md)
- [Training Data](TRAINING_DATA.md)
- [Reproducible production setup](feature_run/ReproducibleSetup.md)
