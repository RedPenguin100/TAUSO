# TAUSO Training Data

This step prepares the antisense-oligonucleotide dataset used to train and test TAUSO's models. It is **not** required to compute features or run off-target analysis — only to train models or run the test suite.

It builds on a completed [data setup](DATA_SETUP.md): canonical-gene assignment uses the genome and the bowtie off-target index, so both must already be in place.

## Pipeline

```bash
# 1. Assign a canonical target gene to every ASO via bulk off-target search
python -m notebooks.data.OligoAI.assign_canonical_gene

# 2. Index the data (assign a stable per-oligo id)
python -m notebooks.utils.data

# 3. Filter and average into the training dataset
python -m notebooks.data.OligoAI.process_data [--cpus N]
```

1. **Canonical genes** — maps each ASO sequence back to a canonical target gene with `find_all_gene_off_targets_BULK`, which is why the genome and bowtie index must be set up first.
2. **Indexing** — assigns each oligo a stable `index_oligo` for downstream joins.
3. **Processing** — standardizes columns, computes structural features, filters out unsupported chemistries, unmapped sequences, and non-standard gapmers, then averages replicate measurements into the dataset the model reads. It also writes diagnostic plots and tables.

On a cluster, run `process_data` with `python -u` (or `PYTHONUNBUFFERED=1`) so logs stream live.
