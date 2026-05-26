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

```bash
# Quick setup (1-2 hours)
tauso setup-all              # Genome, bowtie, omics, raccess
tauso build-cell-context     # Cell-line expression + CAI weights
```

**→ Complete guide:** [Data Setup Guide](documentation/DATA_SETUP.md)

## Functions

Our main function, the `design_asos` function. **NOTE** this function can be slow, with runtime over 30 min. Calculating
the features is not quick in this version of the tool, stay tuned for the article.

# TODO: Improve documentation for the design_asos function.


* Set the `custom_sequence` variable to run an exogenous / mutated gene analysis.
* With `include_details=True`, additional columns will appear with more thorough breakdown. On by default.
* While the API is not settled, many features exist under `tauso.features`, ready to be explored by the user.

## Collaborations

If you want to work with us, or access the future paper, please feel free to contact us at igem.tau.official@gmail.com

