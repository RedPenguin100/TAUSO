# TAUSO Data Setup

TAUSO computes its features from a genome, an off-target alignment index, and several omics datasets. These are downloaded and built locally with the `tauso` CLI rather than shipped with the package. A complete setup is on the order of 10 GB and is dominated by downloading the genome and constructing the bowtie off-target index — there is no shortcut around either if you want off-target analysis.

Install TAUSO first ([Installation](INSTALLATION.md)). Everything below writes to `TAUSO_DATA_DIR` (default `~/.local/share/tauso`); fix that location before you start — see [Choosing a data directory](INSTALLATION.md#choosing-a-data-directory).

## All at once

```bash
tauso setup-all            # genome + bowtie index + omics + raccess
tauso build-cell-context   # per-cell expression, CAI weights, human tGCN
```

`setup-all` is idempotent: datasets already present are verified by hash and skipped unless `--force` is given. On a multi-core machine, pass `-t $(nproc)` to parallelize the bowtie index. Expect the better part of an hour, most of it in indexing.

> **TODO:** `setup-all` does not actually set up everything — it omits `build-cell-context` (above) and the model training data ([Training Data](TRAINING_DATA.md)). Either fold those in, or rename the command so it matches its behaviour.

## One step at a time

### Genome and annotations (required)

```bash
tauso setup-genome [--genome GRCh38]
```

Downloads the FASTA genome and GTF/GFF annotations and builds a gffutils database for fast gene lookups. GRCh38 (human) is supported; GRCm39 (mouse) is a work in progress and not yet fully supported. The FASTA and the annotation database are the largest single items, several GB each.

### Bowtie off-target index (required for off-target analysis)

```bash
tauso setup-bowtie [-t THREADS] [--mem-per-thread MB]
```

Builds the bowtie index used for off-target search. Single-threaded this takes 30–60 minutes and scales with `-t`. Default memory is 800 MB/thread.

### RNA folding tool (required for structure features)

```bash
tauso setup-raccess
```

Compiles `raccess` locally — its license prohibits redistribution, so it cannot be bundled. It needs a C/C++ toolchain and zlib headers; if compilation fails, make sure these are installed:

```bash
# Debian/Ubuntu
sudo apt-get install build-essential zlib1g-dev
# RHEL/CentOS
sudo yum install gcc gcc-c++ zlib-devel
```

### Omics datasets (required for cell-context features)

```bash
tauso setup-omics          # all of the below at once
```

or individually:

```bash
tauso setup-depmap         # DepMap cell-line expression (Public 25Q3), CSV → Parquet
tauso setup-mrna-halflife  # mRNA half-life (TTDB)
tauso setup-tgcn           # tRNA gene copy numbers (GtRNAdb), used for tAI
tauso setup-attract        # ATtRACT RNA-binding-protein motifs
tauso setup-riboseq        # 40S ribosome-profiling track
```

### Cell context (required for per-cell features)

```bash
tauso build-cell-context [CELL_NAMES...] [--reset] [--no-cai]
```

Registers a cohort of cell lines (a default set of 30 if none are named), extracts per-cell expression from DepMap, computes per-cell CAI (Codon Adaptation Index) weights, and sets up the human tGCN data.

```bash
tauso build-cell-context                      # default cohort
tauso build-cell-context HEPG2 HELA A549      # custom cohort (append)
tauso build-cell-context --reset HEPG2 HELA   # replace the cohort
tauso build-cell-context --no-cai             # skip CAI weights
```

Writes `cell_cohort.json`, `processed_expression/`, and `cai_weights.json` to the data directory.

## Verifying

```bash
# Active data directory
python -c "from tauso.data.data import get_data_dir; print(get_data_dir())"

# What has been downloaded
ls -lh "$(python -c 'from tauso.data.data import get_data_dir; print(get_data_dir())')"
```

A complete GRCh38 setup contains, among others, `GRCh38.fa`, `GRCh38.gtf.db`, `GRCh38_bowtie_index/`, the DepMap Parquet files, the mRNA half-life and tGCN tables, `attract/`, the ribo-seq track, `processed_expression/`, `cell_cohort.json`, and `cai_weights.json`.

## Training data

Preparing the model training data is a separate step, needed only to train models or run the test suite. See [Training Data](TRAINING_DATA.md).
