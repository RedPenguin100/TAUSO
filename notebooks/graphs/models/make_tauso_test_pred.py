"""Regenerate the held-out test predictions that Fig 2 reads.

Scores the frozen test split with the released booster (tauso.inference.load_model, md5-pinned)
and writes [index_oligo, tauso] to tauso_test_pred.parquet. The model is trained on train+val, so
it has not seen these rows.

    python make_tauso_test_pred.py
"""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]          # notebooks/graphs/models -> repo root
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
from notebooks.models import common                 # noqa: E402
from tauso.inference import load_model              # noqa: E402

OUT = Path(__file__).resolve().parent / "tauso_test_pred.parquet"


def main():
    df, _ = common.load_dataset()
    _, test = common.split(df)
    model, features = load_model()
    pred = common.predict(model, test, features)
    pd.DataFrame({"index_oligo": test["index_oligo"].to_numpy(), "tauso": pred}).to_parquet(OUT, index=False)
    print(f"wrote {OUT}\n  rows={len(test)}  features={len(features)}")


if __name__ == "__main__":
    main()
