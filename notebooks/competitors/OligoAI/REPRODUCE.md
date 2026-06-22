# Reproducing the TAUSO dataset

Numbered scripts in `notebooks/data/OligoAI/` build the data; raw data and the feature matrix are
pinned on Zenodo. Run from the repo root with the env active and `TAUSO_DATA_DIR` set.

```
1 download → 1.5 split → 2 gene → 2.5 index → 3 filter+average → merge v8 features
  178,624     178,624     178,624   178,624     142,653            train table
```

```bash
python notebooks/data/OligoAI/1_download.py              # raw flank-50 from Zenodo 20794660
python notebooks/data/OligoAI/1_5_assign_split.py        # patent split seed=1  (--verify <ref> optional)
python notebooks/data/OligoAI/2_assign_canonical_gene.py # gene from patent label; NaN if unresolved (kept)
python notebooks/data/OligoAI/2_5_index.py               # add index_oligo (spans the full 178,624 raw)

tauso setup-all && tauso build-cell-context              # genome + omics, once (~10 GB, ~1 h; DATA_SETUP.md)

python notebooks/data/OligoAI/3_process_data.py --cpus $(nproc)   # structure + filter + average → 142,653
```

**Features + merge** — feature matrix = Zenodo 20786540 (`oligo_features_v8.parquet`, auto-fetched or
`tauso setup-features`); `notebooks/models/utility.load_and_validate_final_data` joins it to step 3 on
`index_oligo`.

All filtering lives in step 3 (`notebooks/preprocessing.process_oligo_data`): drops mixmers, missing
gene (incl. the unrecoverable RHO/NAIP), missing cell line, and unmapped. Knobs there: add `"mixmer"`
to the chemistry whitelist → 159,564; `strict_gapmer_patterns=True` → ~131,711. Because `index_oligo`
spans the full raw, the v8 feature parquet must be regenerated whenever the row set changes.

`curate_gene_labels.py` (diagnostic, not in the per-run pipeline) re-derives the gene corrections in
`gene_corrections.py` by aligning ASOs and flagging patents whose `target_gene` disagrees.
