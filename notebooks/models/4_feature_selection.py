"""Step 4 -- feature selection: RMSE-importance ranking + CV K-curve on the chosen objective (train+val only)."""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # repo root, for the notebooks.* imports
from notebooks.models import common
from tauso.data.consts import INHIBITION_PERCENT

VARIANT = "clean_exp"   # default; override with --variant
KGRID = [510, 450, 400, 350, 300, 275, 250, 225, 200, 180, 160, 150, 140,
         120, 100, 90, 80, 70, 60, 50, 40, 35, 30, 25, 20]
STABILITY_KS = [400, 250, 150, 100, 60]   # report ranking stability (fold overlap) at a few K, for whichever you pick
RANK_SEED = 0
PERM_ROWS = 8000


def canon(subset, features):
    """`subset` in the original feature order, so an identical set always yields the same (colsample) model."""
    s = set(subset)
    return [f for f in features if f in s]


def rank_features(train_df, features, units):
    """Permutation-importance ranking by the RMSE increase on the demeaned target (train_df only).
    Per-ASO features are shuffled per row; group-constant features are block-permuted at their level
    (common.permutation_units)."""
    model = common.train(train_df, features, VARIANT, common.PARAMS, common.ROUNDS)
    sample = train_df.sample(n=min(PERM_ROWS, len(train_df)), random_state=RANK_SEED)
    X = np.ascontiguousarray(sample[features].to_numpy(np.float64))
    dy = (sample[INHIBITION_PERCENT]
          - sample.groupby("custom_id")[INHIBITION_PERCENT].transform("mean")).to_numpy()
    rmse = lambda pred: float(np.sqrt(np.mean((pred - dy) ** 2)))
    base = rmse(model.inplace_predict(X))
    rng = np.random.default_rng(RANK_SEED)
    importance = np.empty(len(features))
    for cols, group_col in units:
        saved = X[:, cols].copy()
        common.permute_unit(X, cols, group_col, sample, rng)
        imp = rmse(model.inplace_predict(X)) - base                # permuting an important feature raises RMSE
        X[:, cols] = saved
        for j in cols:
            importance[j] = imp
    return [features[i] for i in np.argsort(-importance)]


def main():
    global VARIANT
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True, choices=list(common.VARIANTS))
    VARIANT = ap.parse_args().variant
    df, features = common.load_dataset()
    trv, _ = common.split(df)                                      # train+val only
    units = common.permutation_units(trv, features)                # per-feature permutation grouping (once)

    fold_rankings, per_fold = [], []
    for tr, va in common.gene_cv_folds(trv):
        ranking = rank_features(tr, features, units)
        fold_rankings.append(ranking)
        scores = {}
        for K in KGRID:
            cols = canon(ranking[:K], features)
            scores[K] = common.metrics_on(common.train(tr, cols, VARIANT, common.PARAMS, common.ROUNDS), va, cols)
        per_fold.append(scores)
    curve = {f"K={K}": common.mean_std([pf[K] for pf in per_fold]) for K in KGRID}

    best_k = max(KGRID, key=lambda K: curve[f"K={K}"]["exp_med"][0])
    best_mean, best_std = curve[f"K={best_k}"]["exp_med"]
    se = best_std / np.sqrt(len(per_fold))
    k_star = min(K for K in KGRID if curve[f"K={K}"]["exp_med"][0] >= best_mean - se)

    final_ranking = rank_features(trv, features, units)            # RMSE ranking over all train+val
    stability = {K: float(np.mean([len(set(fr[:K]) & set(final_ranking[:K])) / K for fr in fold_rankings]))
                 for K in STABILITY_KS}

    lines = [f"Step 4 -- feature selection on [{VARIANT}], nested gene-grouped CV, ranked by RMSE   "
             f"{len(features)} features   train+val={len(trv)}", ""]
    lines += common.metric_table("CV by feature count K  (ranked in-fold by RMSE, scored on held-out fold)", curve) + [""]
    lines += [f"best CV exp_med:  K={best_k}  ({best_mean:.3f} +/- {best_std:.3f}, SE {se:.3f})",
              f"1-SE rule -> K* = {k_star}  (informational; exp_med {curve[f'K={k_star}']['exp_med'][0]:.3f})",
              "",
              "ranking stability (top-K overlap across folds) -- pick any K, take the top-K of feature_ranking.txt:",
              "  " + "   ".join(f"K={K}: {stability[K]:.2f}" for K in STABILITY_KS)]
    common.write_report(lines, f"feature_selection_{VARIANT}.txt")
    common.write_checked(common.RESULTS_DIR / f"feature_ranking_{VARIANT}.txt", "\n".join(final_ranking) + "\n")
    print(f"saved -> feature_ranking_{VARIANT}.txt ({len(final_ranking)} features, ranked most->least important; "
          f"top-K = your K-feature model)")


if __name__ == "__main__":
    main()
