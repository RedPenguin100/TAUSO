"""Objective choice -- compare every model variant (objective x target grouping, see common.VARIANTS) on the
full feature set with the fixed regularized config, by frozen gene-grouped CV on train+val. The winner is
selected by CV exp_med, and feature selection then runs on it.

Run:  python notebooks/models/2_objective_choice.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # repo root, for the notebooks.* imports
from notebooks.models import common

# Regularized config justified in step 1 (default XGBoost overfits these 510 correlated features).
PARAMS = dict(tree_method="hist", device="cuda", max_depth=6, learning_rate=0.05,
              subsample=0.8, colsample_bytree=0.6, min_child_weight=20, reg_lambda=5.0, reg_alpha=1.0)
ROUNDS = 700


def cv_scores(trv, features, variant):
    """Mean +/- std of the metric suite over the frozen gene-grouped CV folds."""
    return common.mean_std([common.metrics_on(common.train(tr, features, variant, PARAMS, ROUNDS), va, features)
                            for tr, va in common.gene_cv_folds(trv)])


def main():
    df, features = common.load_dataset()
    trv, _ = common.split(df)                                      # train+val only

    cv = {v: cv_scores(trv, features, v) for v in common.VARIANTS}
    winner = max(cv, key=lambda v: cv[v]["exp_med"][0])

    lines = [f"Objective choice   {len(features)} features   train+val={len(trv)}   cv=5-fold gene-grouped", ""]
    lines += common.metric_table("CROSS-VALIDATION (train+val, mean +/- std over 5 gene-folds)", cv) + [""]
    lines += [f"winner (by CV exp_med): {winner}  ({cv[winner]['exp_med'][0]:.3f})  <- feature selection runs on this"]
    common.write_report(lines, "objective_choice.txt")


if __name__ == "__main__":
    main()
