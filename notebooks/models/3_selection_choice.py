"""Step 3 -- selection-criterion choice: which feature-importance criterion step 4 ranks by (performance + stability)."""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # repo root, for the notebooks.* imports
from notebooks.models import common
from tauso.data.consts import INHIBITION_PERCENT

VARIANT = "clean_exp"   # default; override with --variant
KGRID = [510, 300, 200, 150, 100, 70, 40]
CRITERIA = ["exp_med", "exp_mean", "wexp_mean", "gxc_med", "rmse"]
HIGHER = {"exp_med": True, "exp_mean": True, "wexp_mean": True, "gxc_med": True, "rmse": False}
PERM_ROWS = 8000
SEED = 0


def canon(subset, features):
    """`subset` in original feature order, so the SET (not order) decides the model (colsample is position-based)."""
    s = set(subset)
    return [f for f in features if f in s]


def criteria(pred, df):
    """All five ranking criteria for `pred` on `df`."""
    y = df[INHIBITION_PERCENT].to_numpy(np.float64)
    dy = y - df.groupby("custom_id")[INHIBITION_PERCENT].transform("mean").to_numpy()
    sp, sz = common.grouped_spearman(pred, y, df["custom_id"].to_numpy(), return_sizes=True)
    gx = common.grouped_spearman(pred, y, df["cohort_id"].to_numpy())
    return {"exp_med": float(np.median(sp)), "exp_mean": float(sp.mean()),
            "wexp_mean": float(np.average(sp, weights=sz)), "gxc_med": float(np.median(gx)),
            "rmse": float(np.sqrt(np.mean((pred - dy) ** 2)))}


def rank_all(train_df, features, units):
    """One model + ONE shared permutation per unit -> importance under EVERY criterion from the same draws,
    so the criteria differ only by metric, never by randomization. Per-ASO features are shuffled per row;
    group-constant features are block-permuted at their level (common.permutation_units). Returns
    {criterion: ordered features}."""
    model = common.train(train_df, features, VARIANT, common.PARAMS, common.ROUNDS)
    sample = train_df.sample(n=min(PERM_ROWS, len(train_df)), random_state=SEED)
    X = np.ascontiguousarray(sample[features].to_numpy(np.float64))
    base = criteria(model.inplace_predict(X), sample)
    imp = {c: np.empty(len(features)) for c in CRITERIA}
    rng = np.random.default_rng(SEED)
    for cols, group_col in units:
        saved = X[:, cols].copy()
        common.permute_unit(X, cols, group_col, sample, rng)     # the SAME permuted unit feeds all criteria
        m = criteria(model.inplace_predict(X), sample)
        X[:, cols] = saved
        for c in CRITERIA:
            val = (base[c] - m[c]) if HIGHER[c] else (m[c] - base[c])
            for j in cols:
                imp[c][j] = val
    return {c: [features[i] for i in np.argsort(-imp[c])] for c in CRITERIA}


def main():
    global VARIANT
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True, choices=list(common.VARIANTS))
    VARIANT = ap.parse_args().variant
    df, features = common.load_dataset()
    trv, _ = common.split(df)                                    # train+val only (no test)
    units = common.permutation_units(trv, features)             # per-feature permutation grouping (once)

    fold_rankings = []
    exp_med = {c: {K: [] for K in KGRID} for c in CRITERIA}
    for tr, va in common.gene_cv_folds(trv):
        rankings = rank_all(tr, features, units)
        fold_rankings.append(rankings)
        for c in CRITERIA:
            for K in KGRID:
                cols = canon(rankings[c][:K], features)
                exp_med[c][K].append(common.metrics_on(common.train(tr, cols, VARIANT, common.PARAMS, common.ROUNDS), va, cols)["exp_med"])
    final = rank_all(trv, features, units)                      # final ranking per criterion (stability reference)

    perf = {c: {K: float(np.mean(exp_med[c][K])) for K in KGRID} for c in CRITERIA}
    stab = {c: {K: float(np.mean([len(set(fr[c][:K]) & set(final[c][:K])) / K for fr in fold_rankings])) for K in KGRID}
            for c in CRITERIA}
    score = {c: float(np.mean([perf[c][K] for K in KGRID])) for c in CRITERIA}      # mean exp_med over the grid
    winner = max(CRITERIA, key=lambda c: score[c])

    head = " ".join(f"K={K:>3}" for K in KGRID)
    lines = [f"Step 3 -- selection-criterion choice on [{VARIANT}], nested gene-grouped CV (no test)", ""]
    lines += ["CV exp_med by ranking criterion  (performance; common yardstick):",
              f"  {'criterion':10s} {head}    mean"]
    for c in CRITERIA:
        lines.append(f"  {c:10s} " + " ".join(f"{perf[c][K]:5.3f}" for K in KGRID) + f"   {score[c]:.3f}")
    lines += ["", "ranking stability -- top-K overlap across folds  (low = a very different model each fold):",
              f"  {'criterion':10s} {head}"]
    for c in CRITERIA:
        lines.append(f"  {c:10s} " + " ".join(f"{stab[c][K]:5.2f}" for K in KGRID))
    lines += ["", f"winner (best mean CV exp_med): {winner}   "
              f"(mean stability {np.mean([stab[winner][K] for K in KGRID]):.2f})  -> step 4 ranks by this"]
    common.write_report(lines, f"selection_choice_{VARIANT}.txt")


if __name__ == "__main__":
    main()
