"""
Per-cell-line and per-modification evaluation across all sweep models + baselines.

Walks `--sweep-dir/*/` for trained TAUSO models, loads data ONCE (with
load_competition=True so baseline scores like oligo_ai_score / PFRED are
columns on every split), then evaluates each model + each baseline on Train /
Val / Test broken down by cell line and by modification, for both custom_id
and gene x cell_line groupings.

Mirrors the analysis in notebooks/models/Evaluation/Evaluation_per_cell_line.ipynb
but is invokable via tauso_run so the heavy data load is amortized across
every model.

Usage:
    tauso_run --cpu=16 --mem=32G \\
      'python -u -m notebooks.models.Evaluation.evaluate_sweep_breakdown \\
         --sweep-dir notebooks/models/Sweep \\
         --out-dir   notebooks/models/Sweep/breakdown'
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

from tauso.data.consts import CANONICAL_GENE, CELL_LINE, INHIBITION, MODIFICATION


COMPETITION = [
    "PFRED_PLS",
    "PFRED_SVM",
    "OW_Overall",
    "OW_Tm",
    "OW_Intra_Oligo",
    "OW_Duplex",
    "sfold_accessibility",
    "miranda_score",
    "miranda_energy",
    "oligo_ai_score",
]


def _large_cohort_indices(df, group_cols, min_size):
    df = df.reset_index(drop=True)
    return [g.index.values for _, g in df.groupby(group_cols) if len(g) >= min_size]


def _metrics(preds, y_true, cohort_idx):
    if len(y_true) == 0:
        return {"Spearman": np.nan, "Top1_Inhib": np.nan, "Top5_Inhib": np.nan,
                "MAE": np.nan, "RMSE": np.nan, "n_cohorts": 0}
    overall_mask = np.isfinite(preds) & np.isfinite(y_true)
    if overall_mask.any():
        mae  = float(mean_absolute_error(y_true[overall_mask], preds[overall_mask]))
        rmse = float(np.sqrt(mean_squared_error(y_true[overall_mask], preds[overall_mask])))
    else:
        mae  = float("nan")
        rmse = float("nan")
    spearmans, top1s, top5s = [], [], []
    for idxs in cohort_idx:
        t, p = y_true[idxs], preds[idxs]
        finite = np.isfinite(t) & np.isfinite(p)
        if finite.sum() < 2:
            continue
        t, p = t[finite], p[finite]
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


def _predict(model_spec, df):
    """Return prediction array for either an XGBoost model or a baseline column."""
    if model_spec["kind"] == "baseline":
        return df[model_spec["col"]].values
    bst   = model_spec["booster"]
    feats = bst.feature_names
    return bst.predict(xgb.DMatrix(df[feats].values, feature_names=feats))


def _discover_models(sweep_dir, model_glob):
    models = []
    for d in sorted(sweep_dir.glob("*/")):
        for m in sorted(d.glob(model_glob)):
            bst = xgb.Booster()
            bst.load_model(str(m))
            models.append({
                "name":       f"{d.name}/{m.stem}",
                "kind":       "xgboost",
                "booster":    bst,
                "n_features": len(bst.feature_names),
            })
    return models


def _baseline_specs(present_in_df):
    specs = []
    for col in COMPETITION:
        if col in present_in_df:
            specs.append({"name": f"Competition ({col})", "kind": "baseline",
                          "col": col, "n_features": 1})
    return specs


def _evaluate(models, splits, breakdown_col, groupings, min_size):
    """
    Returns a long-form DataFrame with columns:
      Model | Split | Breakdown | Grouping | Spearman | Top1_Inhib | Top5_Inhib | MAE | RMSE | n_cohorts | n_features
    `breakdown_col` is None (overall) or the column name used to slice each split.
    """
    rows = []
    for spec in models:
        for split_name, df in splits.items():
            preds  = _predict(spec, df)
            y_true = df[INHIBITION].values

            if breakdown_col is None:
                bins = [(None, np.ones(len(df), dtype=bool))]
            else:
                values = df[breakdown_col].dropna().unique()
                bins = [(v, (df[breakdown_col] == v).values) for v in values]

            for bin_value, mask in bins:
                if not mask.any():
                    continue
                sub = df[mask].reset_index(drop=True)
                p_sub = preds[mask]
                y_sub = y_true[mask]
                for grouping_label, group_cols in groupings.items():
                    idxs = _large_cohort_indices(sub, group_cols, min_size=min_size)
                    m = _metrics(p_sub, y_sub, idxs)
                    rows.append({
                        "Model":      spec["name"],
                        "Split":      split_name,
                        "Breakdown":  "ALL" if bin_value is None else bin_value,
                        "Grouping":   grouping_label,
                        "n_features": spec["n_features"],
                        **m,
                    })
    return pd.DataFrame(rows)


def _print_test_heatmap(df, breakdown_label, grouping_label):
    sub = df[(df["Split"] == "Test") & (df["Grouping"] == grouping_label)]
    if sub.empty:
        print(f"  (no Test rows for grouping={grouping_label})")
        return
    pivot = (sub.pivot_table(index="Model", columns="Breakdown", values="Spearman")
                .round(3))
    # rank models by their median Spearman across breakdowns
    pivot["__median__"] = pivot.median(axis=1)
    pivot = pivot.sort_values("__median__", ascending=False).drop(columns="__median__")
    print(f"\n--- TEST Spearman | breakdown={breakdown_label} | grouping={grouping_label} ---")
    print(pivot.to_string())


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sweep-dir", type=Path, default=Path("notebooks/models/Sweep"),
                        help="Directory containing per-experiment subdirs with trained models")
    parser.add_argument("--extra-model", type=Path, action="append", default=[],
                        help="Path to an additional trained model not under --sweep-dir "
                             "(repeatable). Example: --extra-model "
                             "notebooks/models/SeenOligoModel/Model_Oligo_L2_TrainVal.json")
    parser.add_argument("--model-glob", default="Model_Oligo_*_TrainVal.json",
                        help="Glob for model files inside each experiment dir "
                             "(default: TrainVal -- skips leaked AllData)")
    parser.add_argument("--out-dir", type=Path, default=Path("notebooks/models/Sweep/breakdown"),
                        help="Where to write the CSV outputs")
    parser.add_argument("--min-size-cell-line", type=int, default=20,
                        help="Min cohort size for per-cell-line evaluation (default 20, "
                             "matching the Evaluation_per_cell_line notebook).")
    parser.add_argument("--min-size-modification", type=int, default=1,
                        help="Min cohort size for per-modification evaluation (default 1, "
                             "matching the notebook -- modifications already partition data).")
    args = parser.parse_args()

    # 1. Discover models
    models = _discover_models(args.sweep_dir, args.model_glob)
    for extra in args.extra_model:
        bst = xgb.Booster()
        bst.load_model(str(extra))
        models.append({
            "name":       extra.stem,
            "kind":       "xgboost",
            "booster":    bst,
            "n_features": len(bst.feature_names),
        })

    if not models:
        print(f"No models found under {args.sweep_dir} matching '{args.model_glob}' "
              f"and no --extra-model provided.")
        return

    print(f"Found {len(models)} XGBoost model(s):")
    for m in models:
        print(f"  {m['name']}  ({m['n_features']} features)")

    # 2. Load data ONCE -- with load_competition=True so baselines are present
    print("\nLoading data (with competition baselines)...")
    final_data, _ = load_and_validate_final_data(version="oligo", load_competition=True)
    train_df, val_df, test_df = split_data(
        final_data, [],
        split_col="split", train_val="train", val_val="val", test_val="test",
    )
    print(f"  splits: train={len(train_df)} val={len(val_df)} test={len(test_df)}")

    baselines = _baseline_specs(set(final_data.columns))
    present = [b["col"] for b in baselines]
    missing = [c for c in COMPETITION if c not in present]
    print(f"  baselines present: {present}")
    if missing:
        print(f"  baselines missing (skipped): {missing}")

    all_specs = baselines + models
    splits = {"Train": train_df, "Validation": val_df, "Test": test_df}
    groupings = {
        "custom_id": "custom_id",
        "gene x CL": [CANONICAL_GENE, CELL_LINE],
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # 3. Overall
    print("\n[1/3] Overall (no breakdown)...")
    df_overall = _evaluate(all_specs, splits, None, groupings,
                           min_size=args.min_size_cell_line)
    df_overall.to_csv(args.out_dir / "overall.csv", index=False)

    test_overall = (df_overall[df_overall["Split"] == "Test"]
                    .pivot_table(index="Model", columns="Grouping", values="Spearman")
                    .round(3)
                    .sort_values("gene x CL", ascending=False))
    print("\n--- TEST Spearman | overall ---")
    print(test_overall.to_string())

    # 4. Per cell line
    print("\n[2/3] Per cell line...")
    df_cl = _evaluate(all_specs, splits, CELL_LINE, groupings,
                      min_size=args.min_size_cell_line)
    df_cl.to_csv(args.out_dir / "per_cell_line.csv", index=False)
    _print_test_heatmap(df_cl, "Cell Line", "custom_id")
    _print_test_heatmap(df_cl, "Cell Line", "gene x CL")

    # 5. Per modification
    print("\n[3/3] Per modification...")
    df_mod = _evaluate(all_specs, splits, MODIFICATION, groupings,
                       min_size=args.min_size_modification)
    df_mod.to_csv(args.out_dir / "per_modification.csv", index=False)
    _print_test_heatmap(df_mod, "Modification", "custom_id")
    _print_test_heatmap(df_mod, "Modification", "gene x CL")

    print(f"\nSaved -> {args.out_dir}/overall.csv, per_cell_line.csv, per_modification.csv")


if __name__ == "__main__":
    main()
