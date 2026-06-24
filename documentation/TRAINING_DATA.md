# TAUSO Training Data

This step prepares the antisense-oligonucleotide dataset used to train and test TAUSO's models. It is **not** required to compute features or run off-target analysis — only to train models or run the test suite.

It builds on a completed [data setup](DATA_SETUP.md): the canonical-gene and structure steps use the genome and the bowtie off-target index, so both must already be in place.

## Pipeline

One command runs the whole pipeline — download, split, canonical-gene assignment, indexing, then filter/average:

```bash
python notebooks/data/OligoAI/setup_data.py [--cpus N]
```

It orchestrates the numbered steps in `notebooks/data/OligoAI/`, writing each stage to `raw_data/` (gitignored, reproducible from Zenodo):

1. **`1_download`** — fetch the frozen raw ASO dataset from Zenodo.
2. **`1_5_assign_split`** — reproduce the patent-grouped train/val/test split.
3. **`2_assign_canonical_gene`** — label each ASO with its target gene (patent label plus curated corrections).
4. **`2_5_index`** — assign each oligo a stable `index_oligo` over the full raw set, so downstream joins are order-independent.
5. **`3_process_data`** — compute structure features, filter (unsupported chemistries, unmapped sequences, missing labels), and average replicate measurements into the dataset the model reads.

Step 3 needs the genome (`tauso setup-genome`); pass `--skip-process` to stop after indexing (enough for the feature pipeline, which reads the indexed table). On a cluster, run with `python -u` (or `PYTHONUNBUFFERED=1`) so logs stream live.
