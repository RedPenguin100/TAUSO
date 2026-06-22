"""Step 1.5: assign OligoAI's patent-grouped train/val/test split (seed=1).

Whole patents go to a single split (no patent spans two), reproducing OligoAI's
rinalmo/data/downstream/aso/dataset.py::train_val_test_split. Verified identical to their code.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from notebooks.consts import ORIGINAL_OLIGO_CSV, ORIGINAL_OLIGO_CSV_RAW


def patent_of(custom_id):
    return custom_id.split("/")[-1].split("_table_")[0]


def add_split(df, seed=1, train_ratio=0.8, val_ratio=0.1):
    df = df.dropna(subset=["inhibition_percent"]).reset_index(drop=True)
    patent = df["custom_id"].apply(patent_of)
    np.random.seed(seed)
    shuffled = np.random.permutation(patent.unique())
    n_train = int(train_ratio * len(shuffled))
    n_val = int(val_ratio * len(shuffled))
    train = set(shuffled[:n_train])
    val = set(shuffled[n_train:n_train + n_val])
    df["split"] = np.where(patent.isin(train), "train", np.where(patent.isin(val), "val", "test"))
    return df


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--verify", type=Path, help="Existing split file to assert an exact match against.")
    args = p.parse_args()

    df = add_split(pd.read_csv(ORIGINAL_OLIGO_CSV_RAW))
    print("split counts:", df["split"].value_counts().to_dict())

    if args.verify:
        ref = pd.read_csv(args.verify)
        same = (df["custom_id"].values == ref["custom_id"].values).all() and \
               (df["split"].values == ref["split"].values).all()
        assert same, "reproduced split does not match the reference"
        print("verify: identical to reference")

    df.to_csv(ORIGINAL_OLIGO_CSV, index=False)
    print(f"wrote {len(df):,} rows -> {ORIGINAL_OLIGO_CSV}")


if __name__ == "__main__":
    main()
