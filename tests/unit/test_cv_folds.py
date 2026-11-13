import numpy as np
import pandas as pd
import pytest
from notebooks.data_split.assign_cv_folds import (
    FOLD_COLUMNS,
    INDEX,
    N_FOLDS,
    assign_cv_folds,
    check_cv_folds,
    group_kfold_assignment,
    random_order_assignment,
)

from tauso.data.consts import CANONICAL_GENE_NAME

ROWS_PER_GENE = {"A": 40, "B": 30, "C": 25, "D": 20, "E": 15, "F": 10, "G": 8, "H": 5, "I": 3, "J": 2}


def _table():
    genes = [gene for gene, n in ROWS_PER_GENE.items() for _ in range(n)]
    df = pd.DataFrame({CANONICAL_GENE_NAME: genes, "split": ["train"] * len(genes)})
    df.loc[::3, "split"] = "val"
    test = pd.DataFrame({CANONICAL_GENE_NAME: ["T"] * 4, "split": ["test"] * 4})
    df = pd.concat([df, test], ignore_index=True)
    df[INDEX] = range(1, len(df) + 1)
    return df


@pytest.mark.parametrize(
    "assign",
    [group_kfold_assignment, lambda genes: random_order_assignment(genes, seed=101)],
    ids=["group_kfold", "random_order"],
)
def test_a_gene_is_in_one_fold_and_every_fold_is_used(assign):
    genes = _table()[CANONICAL_GENE_NAME].to_numpy()
    fold = assign(genes)
    assert set(fold) == set(range(N_FOLDS))
    assert pd.Series(fold).groupby(genes).nunique().eq(1).all()


def test_random_order_follows_its_seed():
    genes = _table()[CANONICAL_GENE_NAME].to_numpy()
    assert np.array_equal(random_order_assignment(genes, 101), random_order_assignment(genes, 101))
    assert not np.array_equal(random_order_assignment(genes, 101), random_order_assignment(genes, 202))


def test_test_rows_get_no_fold():
    df = _table()
    folds = assign_cv_folds(df)
    assert list(folds.columns) == [INDEX, *FOLD_COLUMNS]
    assert set(folds[INDEX]) == set(df.loc[df["split"] != "test", INDEX])
    check_cv_folds(folds, df)


def test_check_rejects_a_gene_split_across_folds():
    df = _table()
    folds = assign_cv_folds(df)
    first_a = df.index[df[CANONICAL_GENE_NAME] == "A"][0]
    row = folds[INDEX] == df.loc[first_a, INDEX]
    folds.loc[row, "fold1"] = (folds.loc[row, "fold1"] + 1) % N_FOLDS
    with pytest.raises(AssertionError, match="fold1"):
        check_cv_folds(folds, df)
