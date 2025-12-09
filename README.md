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

## Installation for python users

The best way to use the library is with a separate conda environment

* `conda env create -f environment.yml`
* `conda activate tauso`
* `pip install -e .`

## Commands and setup

The tool allows you to install genomes easily and index them. Every command has a `--help` option.

`tauso install-raccess` - Install a folding prediction tool named raccess. Distribution of the code is prohibited, so
the user has to install manually. **This command is necessary to run the tool**

`tauso setup-genome` - download the genome and annotation files. Currently, GRCh38 and GRCm39 are supported. Can take
5-10+ minutes. **This command is necessary to run the tool**

`tauso setup-bowtie` - after running `setup-genome`, you may run this command to index the genome. May take 30-60+
minutes.

`tauso run-off-target` - After running both `setup-genome` and `setup-bowtie`, you may analyze off-targets and limit the
mismatches, with useful annotation about each off-target (intergenic, exonic, intronic regions).

NOTE: In the future we plan to write an article about better off-target analysis, this will be posted in a separate
paper.

## Functions

Our main function, the `design_asos` function. **NOTE** this function can be slow, with runtime over 30 min. Calculating
the features is not quick in this version of the tool, stay tuned for the article.

```python

from tauso.old_model_generation.main import design_asos

aso_results, off_target_analysis = design_asos(gene_name='DDX11L1', organism='human', top_k=5, run_off_target=True)

```
* Set the `custom_sequence` variable to run an exogenous / mutated gene analysis.
* With `include_details=True`, additional columns will appear with more thorough breakdown. On by default.
* While the API is not settled, many features exist under `tauso.features`, ready to be explored by the user.

## Collaborations

If you want to work with us, or access the future paper, please feel free to contact us at igem.tau.official@gmail.com

