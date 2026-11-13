"""Assign the gene-grouped cross-validation folds used to compare model configurations.

The train+val rows of the OligoAI split are divided into five folds three different ways, always by
canonical gene, so that no gene is on both the training side and the validation side of a fold. A
configuration is scored by averaging over the 3 x 5 = 15 fits (partition p, fold f: validate on the rows
with fold<p> == f, train on the rest).

    fold0  scikit-learn GroupKFold(5): genes by decreasing size, each to the fold with the fewest rows so far
    fold1  genes in a random order (numpy default_rng(101)), each to the fold with the fewest rows so far
    fold2  the same with default_rng(202)

GroupKFold breaks ties between genes of equal size with np.argsort, so fold0 depends on the numpy version;
cv_folds.csv.gz next to this file is the fold assignment that the model selection reads.

    python notebooks/data_split/assign_cv_folds.py [--table CSV] [--out FILE]   # assign afresh
    python notebooks/data_split/assign_cv_folds.py --check FILE [--table CSV]   # validate a fold file
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from notebooks.consts import OLIGO_AI_DATA_PATH, OLIGO_CSV_PROCESSED_AVERAGED
from tauso.data.consts import CANONICAL_GENE_NAME

INDEX = "index_oligo"
N_FOLDS = 5
RANDOM_ORDER_SEEDS = (101, 202)
FOLD_COLUMNS = [f"fold{p}" for p in range(1 + len(RANDOM_ORDER_SEEDS))]
FROZEN_FOLDS = Path(__file__).resolve().parent / "cv_folds.csv.gz"


def train_val_rows(df):
    return df[df["split"].isin(["train", "val"])].reset_index(drop=True)


def group_kfold_assignment(genes, n_folds=N_FOLDS):
    """Fold of every row from scikit-learn GroupKFold grouped by gene."""
    fold = np.empty(len(genes), dtype=int)
    for k, (_, validation) in enumerate(GroupKFold(n_folds).split(np.zeros(len(genes)), groups=genes)):
        fold[validation] = k
    return fold


def random_order_assignment(genes, seed, n_folds=N_FOLDS):
    """Fold of every row: genes taken in a random order, each to the fold with the fewest rows so far."""
    unique, rows_per_gene = np.unique(genes, return_counts=True)
    size = dict(zip(unique, rows_per_gene))
    rng = np.random.default_rng(seed)
    load = np.zeros(n_folds)
    fold_of_gene = {}
    for gene in rng.permutation(unique):
        k = int(np.argmin(load))
        fold_of_gene[gene] = k
        load[k] += size[gene]
    return np.array([fold_of_gene[gene] for gene in genes])


def assign_cv_folds(df):
    """index_oligo and fold0, fold1, fold2 for every train+val row of `df`."""
    trv = train_val_rows(df)
    assert trv[CANONICAL_GENE_NAME].notna().all(), "every train+val row needs a canonical gene"
    genes = trv[CANONICAL_GENE_NAME].to_numpy()
    folds = pd.DataFrame({INDEX: trv[INDEX].to_numpy(), "fold0": group_kfold_assignment(genes)})
    for column, seed in zip(FOLD_COLUMNS[1:], RANDOM_ORDER_SEEDS):
        folds[column] = random_order_assignment(genes, seed)
    return folds


def check_cv_folds(folds, df):
    """Assert that `folds` covers exactly the train+val rows of `df` and that no gene spans two folds."""
    trv = train_val_rows(df)
    assert list(folds.columns) == [INDEX, *FOLD_COLUMNS], f"unexpected columns {list(folds.columns)}"
    assert folds[INDEX].is_unique, f"{INDEX} repeats in the fold file"
    assert set(folds[INDEX]) == set(trv[INDEX]), "the fold file and the table do not have the same train+val rows"
    merged = trv[[INDEX, CANONICAL_GENE_NAME]].merge(folds, on=INDEX)
    for column in FOLD_COLUMNS:
        assert set(merged[column]) == set(range(N_FOLDS)), f"{column} does not use folds 0-{N_FOLDS - 1}"
        spanning = merged.groupby(CANONICAL_GENE_NAME)[column].nunique().gt(1).sum()
        assert spanning == 0, f"{column}: {spanning} genes are in more than one fold"


def fold_sizes(folds):
    return pd.DataFrame({column: folds[column].value_counts().sort_index() for column in FOLD_COLUMNS})


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--table", type=Path, default=OLIGO_CSV_PROCESSED_AVERAGED, help="OligoAI processed table.")
    p.add_argument("--out", type=Path, default=OLIGO_AI_DATA_PATH / "cv_folds.csv.gz", help="Where to write folds.")
    p.add_argument(
        "--check",
        type=Path,
        nargs="?",
        const=FROZEN_FOLDS,
        help="Validate this fold file (default: the one in this folder) instead of assigning.",
    )
    args = p.parse_args()

    df = pd.read_csv(args.table, usecols=[INDEX, "split", CANONICAL_GENE_NAME])
    if args.check:
        folds = pd.read_csv(args.check)
        check_cv_folds(folds, df)
        print(f"{args.check}: ok, {len(folds):,} train+val rows, no gene in more than one fold")
    else:
        folds = assign_cv_folds(df)
        check_cv_folds(folds, df)
        folds.to_csv(args.out, index=False, compression={"method": "gzip", "mtime": 0})
        print(f"assigned {len(folds):,} train+val rows -> {args.out}")
    print("rows per fold:")
    print(fold_sizes(folds).to_string())


if __name__ == "__main__":
    main()
