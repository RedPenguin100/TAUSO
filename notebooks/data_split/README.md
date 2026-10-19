# Cross-validation folds

The folds used to compare model configurations (hyperparameter search) on the OligoAI train+val rows.
Every fold is grouped by canonical gene, so no gene is on both the training and the validation side.

`assign_cv_folds.py` holds the assignment rules. `cv_folds.csv.gz` (`index_oligo`, `fold0`, `fold1`, `fold2`)
is the assignment the model selection read.

| column | rule |
|---|---|
| `fold0` | scikit-learn `GroupKFold(5)` by canonical gene |
| `fold1` | genes in a random order (`numpy.random.default_rng(101)`), each gene to the fold with the fewest rows so far |
| `fold2` | the same with `default_rng(202)` |

A configuration is scored by 3 partitions x 5 folds = 15 fits: for partition `p` and fold `f`, validate on the
rows with `fold<p> == f` and train on the other rows, then average each metric over the 15 fits.

```bash
python notebooks/data_split/assign_cv_folds.py --check      # validate cv_folds.csv.gz against the processed table
python notebooks/data_split/assign_cv_folds.py              # assign afresh -> raw_data/cv_folds.csv.gz (gitignored)
```

## Which folds are in `cv_folds.csv.gz`

- They were assigned on an earlier version of the processed table (121,997 train+val rows, 426 more than now).
  The file keeps the folds of the 121,571 train+val rows of the current table; the 426 other rows are never
  read, and no gene is in two folds under the current gene labels.
- Assigning afresh on the current table gives a different, equally valid partition (it moves 32,916 / 73,941 /
  40,523 of the 121,571 rows in `fold0` / `fold1` / `fold2`). The file is the record of the folds that were used.
- The functions in `assign_cv_folds.py` give the earlier folds back exactly (0 of 121,997 rows differ) when run
  on the earlier tables, with numpy 2.2.6 and 2.5.4 and scikit-learn 1.7.2.
- `GroupKFold` breaks ties between genes of equal size with `np.argsort`, whose tie order changed in numpy 2:
  with numpy 1.26 the same code gives a different `fold0`.
- `gene_cv_folds` in `notebooks/models/common.py` is a `GroupKFold(5)` on the current rows. It assigns genes to
  folds exactly as a fresh `fold0` does, which is not the `fold0` in `cv_folds.csv.gz`.
