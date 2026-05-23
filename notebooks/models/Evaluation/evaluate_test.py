"""
Evaluate a trained XGBoost ASO model on train/val/test splits.

Reports Spearman, Top1%/Top5% median inhibition, MAE, RMSE — for each split,
grouped both by custom_id (within-batch, easy/inflated) and gene x cell-line
(cross-cohort, the honest generalization metric).

Usage:
    python -u -m notebooks.models.Evaluation.evaluate_test \
        --model notebooks/models/SeenOligoModel/Model_Oligo_L2_TrainVal.json

    python -u -m notebooks.models.Evaluation.evaluate_test \
        --model notebooks/models/SeenOligoModel/Model_Oligo_L2_TrainVal.json \
        --min-eval-size 20 \
        --save-json /tmp/test_eval.json

Run inside the tauso_claude env (mamba run -n tauso_claude ...).
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error

from notebooks.models.SeenOligoModel.base_model import split_data
from notebooks.models.utility import load_and_validate_final_data

from tauso.data.consts import CANONICAL_GENE, CELL_LINE, INHIBITION


def _large_cohort_indices(df, group_cols, min_size):
    df = df.reset_index(drop=True)
    return [g.index.values for _, g in df.groupby(group_cols) if len(g) >= min_size]


def _metrics(preds, y_true, cohort_idx):
    mae = float(mean_absolute_error(y_true, preds))
    rmse = float(np.sqrt(mean_squared_error(y_true, preds)))
    spearmans, top1s, top5s = [], [], []
    for idxs in cohort_idx:
        t, p = y_true[idxs], preds[idxs]
        corr, _ = spearmanr(t, p)
        if not np.isnan(corr):
            spearmans.append(corr)
        n = len(t)
        k1, k5 = max(1, int(n * 0.01)), max(1, int(n * 0.05))
        top1s.append(np.median(t[np.argpartition(p, -k1)[-k1:]]))
        top5s.append(np.median(t[np.argpartition(p, -k5)[-k5:]]))
    return {
        "Spearman":   float(np.nanmedian(spearmans)) if spearmans else float("nan"),
        "Top1_Inhib": float(np.nanmedian(top1s))     if top1s     else float("nan"),
        "Top5_Inhib": float(np.nanmedian(top5s))     if top5s     else float("nan"),
        "MAE":        mae,
        "RMSE":       rmse,
        "n_cohorts":  len(spearmans),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True, type=Path,
                        help="Path to the saved XGBoost JSON model")
    parser.add_argument("--min-eval-size", type=int, default=50,
                        help="Min cohort size to include in median (default: 50). "
                             "Smaller values include more cohorts but noisier estimates.")
    parser.add_argument("--save-json", type=Path, default=None,
                        help="Optional path to dump the metrics rows as JSON")
    args = parser.parse_args()

    is_alldata = "AllData" in args.model.name
    print(f"Loading model: {args.model}")
    bst = xgb.Booster()
    bst.load_model(str(args.model))
    feats = bst.feature_names
    print(f"  features: {len(feats)}")
    if is_alldata:
        print("  WARNING: '_AllData' in filename — train/val/test all seen during training. Test is leaked.")

    print("Loading data...")
    final_data, _ = load_and_validate_final_data(version="oligo")

    train_df, val_df, test_df = split_data(
        final_data, feats,
        split_col="split", train_val="train", val_val="val", test_val="test",
    )
    print(f"  splits: train={len(train_df)} val={len(val_df)} test={len(test_df)}")

    groupings = {
        "custom_id (within-batch)":     "custom_id",
        "gene x cell_line (cross)":    [CANONICAL_GENE, CELL_LINE],
    }

    rows = []
    for split_name, df in [("Train", train_df), ("Validation", val_df), ("Test", test_df)]:
        preds = bst.predict(xgb.DMatrix(df[feats].values, feature_names=feats))
        y_true = df[INHIBITION].values
        for grouping_label, group_cols in groupings.items():
            idxs = _large_cohort_indices(df, group_cols, min_size=args.min_eval_size)
            m = _metrics(preds, y_true, idxs)
            rows.append({"Split": split_name, "Grouping": grouping_label, **m})

    df_results = pd.DataFrame(rows)

    print(f"\n=== Metrics (min_eval_size={args.min_eval_size}) ===")
    print(df_results.to_string(index=False,
        formatters={
            "Spearman":   lambda x: f"{x:.4f}",
            "Top1_Inhib": lambda x: f"{x:.2f}",
            "Top5_Inhib": lambda x: f"{x:.2f}",
            "MAE":        lambda x: f"{x:.3f}",
            "RMSE":       lambda x: f"{x:.3f}",
        }))

    if args.save_json:
        args.save_json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.save_json, "w") as f:
            json.dump(rows, f, indent=2)
        print(f"\nSaved -> {args.save_json}")


if __name__ == "__main__":
    main()
