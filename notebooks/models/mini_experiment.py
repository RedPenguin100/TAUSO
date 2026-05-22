"""
Mini experiment: compare L2, Huber, and rank:pairwise on ~150 ASOs.

Loads a small stratified sample (a few cohorts, fixed row count), trains
three models with fixed hyperparams (no Optuna), and reports per-cohort
Spearman + MAE side by side.

Usage:
    python -m notebooks.models.mini_experiment
    python -m notebooks.models.mini_experiment --n-cohorts 5 --n-per-cohort 30 --cpus 4
"""
import argparse
import json
import logging
import sys

import numpy as np
import pandas as pd
import xgboost as xgb
from scipy.stats import spearmanr

from notebooks.models.utility import load_and_validate_final_data
from tauso.data.consts import CANONICAL_GENE, CELL_LINE, INHIBITION

logger = logging.getLogger(__name__)

# Fixed hyperparams — no Optuna, just reasonable defaults for a small dataset
BASE_PARAMS = {
    "tree_method": "hist",
    "device": "cpu",
    "max_depth": 4,
    "learning_rate": 0.05,
    "subsample": 0.8,
    "colsample_bytree": 0.5,
    "min_child_weight": 5,
    "reg_alpha": 1.0,
    "reg_lambda": 1.0,
    "seed": 1,
}
NUM_ROUNDS = 300
EARLY_STOP = 30


def _cohort_spearman(y_true, preds, groups):
    scores = []
    for idxs in groups:
        if len(idxs) < 3:
            continue
        r, _ = spearmanr(y_true[idxs], preds[idxs])
        if not np.isnan(r):
            scores.append(r)
    return float(np.median(scores)) if scores else float("nan"), len(scores)


def _pick_cohorts(df, n_cohorts, n_per_cohort, split_col="split"):
    """
    Pick n_cohorts that have enough train AND val rows.
    Return a DataFrame subset with at most n_per_cohort train rows per cohort
    (all val rows kept).
    """
    cohort_key = [CANONICAL_GENE, CELL_LINE]
    train_df = df[df[split_col] == "train"]
    val_df   = df[df[split_col] == "val"]

    train_sizes = train_df.groupby(cohort_key).size()
    val_sizes   = val_df.groupby(cohort_key).size()

    eligible = (
        train_sizes[train_sizes >= n_per_cohort]
        .index.intersection(val_sizes[val_sizes >= 5].index)
    )
    if len(eligible) < n_cohorts:
        raise ValueError(f"Only {len(eligible)} cohorts qualify (need {n_cohorts}). "
                         "Lower --n-per-cohort or --n-cohorts.")

    # Pick cohorts spread across different genes
    chosen = list(eligible[:n_cohorts])
    mask = df[cohort_key].apply(tuple, axis=1).isin([tuple(c) for c in chosen])
    subset = df[mask].copy()

    # Cap train rows per cohort
    train_rows = []
    for cohort in chosen:
        cdf = subset[(subset[CANONICAL_GENE] == cohort[0]) &
                     (subset[CELL_LINE] == cohort[1]) &
                     (subset[split_col] == "train")]
        train_rows.append(cdf.sample(min(n_per_cohort, len(cdf)), random_state=1))
    val_rows = subset[subset[split_col] == "val"]

    return pd.concat(train_rows + [val_rows]).reset_index(drop=True), chosen


def _cohort_groups(df, features):
    """Return (X, y, group_sizes, cohort_indices) sorted by cohort."""
    cohort_key = [CANONICAL_GENE, CELL_LINE]
    df = df.sort_values(cohort_key).reset_index(drop=True)
    groups, cohort_indices = [], []
    for _, g in df.groupby(cohort_key, sort=False):
        groups.append(len(g))
        cohort_indices.append(g.index.values)
    X = df[features].values.astype(np.float32)
    y = df[INHIBITION].values.astype(np.float32)
    return X, y, groups, cohort_indices


def _train_regression(objective, X_tr, y_tr, X_vl, y_vl, nthread, extra_params=None):
    params = {**BASE_PARAMS, "objective": objective, "nthread": nthread}
    if extra_params:
        params.update(extra_params)
    dtrain = xgb.DMatrix(X_tr, label=y_tr)
    dval   = xgb.DMatrix(X_vl, label=y_vl)
    bst = xgb.train(
        params, dtrain, num_boost_round=NUM_ROUNDS,
        evals=[(dval, "val")], early_stopping_rounds=EARLY_STOP,
        verbose_eval=False,
    )
    return bst.predict(xgb.DMatrix(X_vl))


def _train_pairwise(X_tr, y_tr, groups_tr, X_vl, y_vl, nthread):
    params = {**BASE_PARAMS, "objective": "rank:pairwise", "nthread": nthread}
    # For ranking, DMatrix needs group sizes
    dtrain = xgb.DMatrix(X_tr, label=y_tr)
    dtrain.set_group(groups_tr)
    # Val: use standard DMatrix (no group needed for prediction)
    dval_rank = xgb.DMatrix(X_vl, label=y_vl)
    bst = xgb.train(
        params, dtrain, num_boost_round=NUM_ROUNDS,
        verbose_eval=False,  # rank:pairwise has no standard eval metric with DMatrix
    )
    return bst.predict(xgb.DMatrix(X_vl))


def run(args):
    logger.info("Loading data (this may take ~30s)...")
    final_data, features = load_and_validate_final_data(version="oligo")
    logger.info("Loaded %d rows, %d features.", len(final_data), len(features))

    logger.info("Sampling cohorts: %d cohorts × up to %d train rows each...",
                args.n_cohorts, args.n_per_cohort)
    subset, chosen = _pick_cohorts(final_data, args.n_cohorts, args.n_per_cohort)
    logger.info("Chosen cohorts: %s", [(c[0], c[1]) for c in chosen])

    train_df = subset[subset["split"] == "train"]
    val_df   = subset[subset["split"] == "val"]
    logger.info("Train rows: %d | Val rows: %d", len(train_df), len(val_df))

    X_tr = train_df[features].values.astype(np.float32)
    y_tr = train_df[INHIBITION].values.astype(np.float32)
    X_vl = val_df[features].values.astype(np.float32)
    y_vl = val_df[INHIBITION].values.astype(np.float32)

    # Cohort indices for val evaluation
    _, _, groups_tr, _ = _cohort_groups(train_df, features)
    _, _, _, val_cohort_idx = _cohort_groups(val_df, features)

    def report(name, preds):
        spear, n_cohorts = _cohort_spearman(y_vl, preds, val_cohort_idx)
        mae = float(np.mean(np.abs(y_vl - preds)))
        logger.info("%-20s | val_spearman=%.4f (over %d cohorts) | MAE=%.3f",
                    name, spear, n_cohorts, mae)
        return {"spearman": spear, "mae": mae, "n_cohorts": n_cohorts}

    results = {}

    logger.info("\n--- L2 (reg:squarederror) ---")
    preds_l2 = _train_regression("reg:squarederror", X_tr, y_tr, X_vl, y_vl, args.cpus)
    results["L2"] = report("L2", preds_l2)

    logger.info("\n--- Huber (reg:pseudohubererror, slope=2.0) ---")
    preds_huber = _train_regression("reg:pseudohubererror", X_tr, y_tr, X_vl, y_vl,
                                    args.cpus, extra_params={"huber_slope": 2.0})
    results["Huber"] = report("Huber (slope=2.0)", preds_huber)

    logger.info("\n--- Pairwise ranking (rank:pairwise) ---")
    preds_pw = _train_pairwise(X_tr, y_tr, groups_tr, X_vl, y_vl, args.cpus)
    results["Pairwise"] = report("Pairwise", preds_pw)

    logger.info("\n=== SUMMARY ===")
    for name, r in results.items():
        logger.info("%-20s | spearman=%.4f | MAE=%.3f", name, r["spearman"], r["mae"])

    out = args.output
    if out:
        import json
        with open(out, "w") as f:
            json.dump(results, f, indent=2)
        logger.info("Results saved -> %s", out)

    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n-cohorts",    type=int, default=5,
                        help="Number of gene×cell-line cohorts to sample (default: 5)")
    parser.add_argument("--n-per-cohort", type=int, default=30,
                        help="Max train rows per cohort (default: 30)")
    parser.add_argument("--cpus",         type=int, default=4,
                        help="XGBoost nthread (default: 4)")
    parser.add_argument("--output",       default=None,
                        help="Optional path to save results JSON")
    parser.add_argument("--log-level",    default="INFO",
                        choices=["DEBUG", "INFO", "WARNING"])
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s | %(message)s",
        stream=sys.stdout, force=True,
    )
    run(args)


if __name__ == "__main__":
    main()
