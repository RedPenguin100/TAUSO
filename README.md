This work is NOT published as-of yet, and is in progress. Please stay tuned for the published version of our software.



Shield: [![CC BY-NC 4.0][cc-by-nc-shield]][cc-by-nc]

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/

[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png

[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg 

# TAUSO

Feature extraction and analysis utilities for antisense oligonucleotides (ASO) design, built for the TAU-Israel 2025
iGEM project.

## Features

* MOE (20-mer) and LNA (16-mer) candidate ranking pipeline with feature breakdown and off-target analysis.
* Features include GC content, hybridization, RNA accessibility, folding, and toxicity heuristics.
* Tool to download genomic sequences, annotations, and index them with `bowtie`.

## Getting Started

**→ New users: See [Quick Start Guide](documentation/QUICKSTART.md) for step-by-step instructions.**

## Installation

```bash
# One command - installs everything
mamba env create -f environment-dev.yml
mamba activate tauso
```

**→ Complete guide:** [Installation Guide](documentation/INSTALLATION.md)

## Data Setup

Features are computed against a locally built genome, off-target index, and omics datasets — on the order of 10 GB,
most of it the genome and the bowtie off-target index. The build is dominated by downloading the genome and indexing
it; budget the better part of an hour.

```bash
tauso setup-all              # genome, bowtie off-target index, omics, raccess
tauso build-cell-context     # per-cell expression + CAI weights + tGCN
```

**→ Complete guide:** [Data Setup Guide](documentation/DATA_SETUP.md)

## Usage

The pipeline ranks MOE (20-mer) and LNA (16-mer) candidates against a target, with a per-feature breakdown and
off-target analysis. Feature computation is the bottleneck in this version — a full target can take tens of minutes —
and performance work is ongoing.

The public API is not yet settled. Feature implementations live under `tauso.features` and can be explored directly.
Pass a `custom_sequence` to analyze an exogenous or mutated target. A documented, stable entry point will accompany
the article.

## Collaborations

If you want to work with us, or access the future paper, please feel free to contact us at igem.tau.official@gmail.com

