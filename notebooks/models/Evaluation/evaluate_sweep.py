"""
Evaluate ALL finished models in a sweep directory.

Walks `--sweep-dir/*/` looking for `Model_Oligo_*_TrainVal.json` (configurable),
loads the data ONCE, then scores every found model on train/val/test for both
custom_id and gene x cell_line groupings. Prints a compact test-Spearman
comparison plus a full per-split table. Optionally dumps a CSV.

Designed to be invoked once via SLURM (tauso_run) so the slow data load is
amortized across all experiments.

Usage:
    python -u -m notebooks.models.Evaluation.evaluate_sweep \
        --sweep-dir notebooks/models/Sweep \
        --min-eval-size 50 \
        --save-csv notebooks/models/Sweep/results.csv
"""
import argparse
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
    mae  = float(mean_absolute_error(y_true, preds))
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
    parser.add_argument("--sweep-dir", type=Path, default=Path("notebooks/models/Sweep"),
                        help="Directory containing per-experiment subdirs with trained models")
    parser.add_argument("--model-glob", default="Model_Oligo_*_TrainVal.json",
                        help="Glob to find model files within each experiment dir "
                             "(default: TrainVal only -- skips leaked AllData)")
    parser.add_argument("--min-eval-size", type=int, default=50,
                        help="Minimum cohort size for median (default: 50)")
    parser.add_argument("--save-csv", type=Path, default=None,
                        help="Optional path to write the full results as CSV")
    args = parser.parse_args()

    # 1. Discover finished experiments
    models = []
    for d in sorted(args.sweep_dir.glob("*/")):
        for m in sorted(d.glob(args.model_glob)):
            models.append((d.name, m))

    if not models:
        print(f"No models found under {args.sweep_dir} matching '{args.model_glob}'")
        print("Nothing has finished yet, or the sweep directory layout is different.")
        return

    print(f"Found {len(models)} finished model(s):")
    for exp, m in models:
        print(f"  {exp}/{m.name}")
    print()

    # 2. Load data ONCE
    print("Loading data (one-time cost)...")
    final_data, _ = load_and_validate_final_data(version="oligo")
    train_df, val_df, test_df = split_data(
        final_data, [],   # features arg unused by split_data
        split_col="split", train_val="train", val_val="val", test_val="test",
    )
    print(f"  splits: train={len(train_df)} val={len(val_df)} test={len(test_df)}")

    groupings = {
        "custom_id": "custom_id",
        "gene x CL": [CANONICAL_GENE, CELL_LINE],
    }

    # 3. Pre-compute cohort indices per split per grouping (model-independent)
    indices = {}
    for split_name, df in [("Train", train_df), ("Validation", val_df), ("Test", test_df)]:
        for grouping_label, group_cols in groupings.items():
            indices[(split_name, grouping_label)] = _large_cohort_indices(
                df, group_cols, min_size=args.min_eval_size,
            )

    # 4. Evaluate every model
    rows = []
    for exp_name, model_path in models:
        print(f"\n=== {exp_name} :: {model_path.name} ===")
        bst = xgb.Booster()
        bst.load_model(str(model_path))
        feats = bst.feature_names
        print(f"  features: {len(feats)}")

        for split_name, df in [("Train", train_df), ("Validation", val_df), ("Test", test_df)]:
            preds = bst.predict(xgb.DMatrix(df[feats].values, feature_names=feats))
            y_true = df[INHIBITION].values
            for grouping_label in groupings:
                m = _metrics(preds, y_true, indices[(split_name, grouping_label)])
                rows.append({
                    "Experiment": exp_name,
                    "Model":      model_path.name,
                    "Split":      split_name,
                    "Grouping":   grouping_label,
                    "n_features": len(feats),
                    **m,
                })

    df_results = pd.DataFrame(rows)

    # 5. Compact test-only Spearman comparison
    test_only = df_results[df_results["Split"] == "Test"]
    test_pivot = (test_only
                  .pivot_table(index=["Experiment", "n_features"],
                               columns="Grouping",
                               values="Spearman")
                  .round(4)
                  .sort_values("gene x CL", ascending=False))

    print("\n" + "=" * 60)
    print("TEST SPEARMAN (ranked by gene x CL)")
    print("=" * 60)
    print(test_pivot.to_string())

    # 6. Full per-split table
    print("\n" + "=" * 60)
    print("FULL RESULTS")
    print("=" * 60)
    print(df_results[["Experiment", "Split", "Grouping", "Spearman",
                      "Top1_Inhib", "MAE", "n_features", "n_cohorts"]]
          .to_string(index=False, formatters={
              "Spearman":   lambda x: f"{x:.4f}",
              "Top1_Inhib": lambda x: f"{x:.2f}",
              "MAE":        lambda x: f"{x:.3f}",
          }))

    if args.save_csv:
        args.save_csv.parent.mkdir(parents=True, exist_ok=True)
        df_results.to_csv(args.save_csv, index=False)
        print(f"\nSaved -> {args.save_csv}")


if __name__ == "__main__":
    main()
