"""Comprehensive LambdaRank vs L2 probe on the relaxed dataset with masked chem features.

Same v3 wide cache (685 cols) + 17 new chem features (with the NaN-masking baked in
on this branch). Relaxed-mode processed_averaged CSV: 142,653 averaged rows.
COMPETITION outputs are excluded from features.

For each variant we report mean + median Spearman, NDCG@5, NDCG@10, P@5, P@10,
Top5_Inhib, Top10_Inhib, MAE -- on BOTH custom_id and gene x cell-line cohorts.

Variants:

  baseline:
    L2                          reg:squarederror, no group, ROUNDS=700

  pairwise:
    pair_custid                 rank:pairwise, group=custom_id, ROUNDS=700
    pair_gene_cell              rank:pairwise, group=gene+cell, ROUNDS=700
    pair_custid_long            rank:pairwise, group=custom_id, ROUNDS=1400

  ndcg (lambdarank with NDCG gain):
    ndcg_custid_default         rank:ndcg, group=custom_id, default pair method
    ndcg_custid_mean_8          rank:ndcg, group=custom_id, mean sampling, 8 pairs/sample
    ndcg_custid_topk_4          rank:ndcg, group=custom_id, topk sampling, 4 pairs/sample
    ndcg_gene_cell_default      rank:ndcg, group=gene+cell, default pair method
    ndcg_gene_cell_mean_8       rank:ndcg, group=gene+cell, mean sampling, 8 pairs/sample
"""

import argparse
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, ndcg_score

from notebooks.consts import OLIGO_CSV_PROCESSED_AVERAGED
from notebooks.features.feature_extraction import COMPETITION
from notebooks.models._chem_overhaul_ab import ADDED_IN_NEW, DROPPED_IN_NEW, add_new_features
from tauso.data.consts import CANONICAL_GENE, CELL_LINE, INHIBITION, VOLUME
from tauso.populate.feature_cache import ensure_cache, index_col as cache_index_col

logger = logging.getLogger(__name__)

SANE_BASE = {
    "tree_method": "hist", "device": "cuda",
    "max_depth": 6, "learning_rate": 0.05,
    "subsample": 0.8, "colsample_bytree": 0.6, "min_child_weight": 20,
    "reg_lambda": 5.0, "reg_alpha": 1.0,
}
ROUNDS_DEFAULT = 700
SEEDS = (1, 2, 3)


def _group_array(df, group_choice):
    if group_choice == "custom_id":
        keys = ["custom_id"]
    elif group_choice == "gene_cell":
        keys = [CANONICAL_GENE, CELL_LINE]
    else:
        raise ValueError(group_choice)
    sorted_df = df.sort_values(keys, kind="stable").reset_index(drop=True)
    sizes = sorted_df.groupby(keys, sort=False).size().to_numpy()
    assert sizes.sum() == len(sorted_df)
    return sorted_df, sizes


def _discretise_for_ndcg(y, n_levels=32):
    order = np.argsort(np.argsort(y, kind="stable"), kind="stable")
    return ((order / max(len(order) - 1, 1)) * (n_levels - 1)).round().astype(np.int32)


def _train(trv, feats, seed, *, objective, group_choice=None,
           rounds=ROUNDS_DEFAULT, lambdarank_pair_method=None, lambdarank_num_pair=None):
    params = {**SANE_BASE, "objective": objective, "seed": seed}

    if objective == "reg:squarederror":
        y = trv[INHIBITION].to_numpy(np.float64)
        X = trv[feats].to_numpy(np.float32)
        d = xgb.QuantileDMatrix(X, label=y, feature_names=feats)
        return xgb.train(params, d, num_boost_round=rounds)

    sorted_trv, groups = _group_array(trv, group_choice)
    y_raw = sorted_trv[INHIBITION].to_numpy(np.float64)
    X = sorted_trv[feats].to_numpy(np.float32)
    if objective == "rank:ndcg":
        y = _discretise_for_ndcg(y_raw, n_levels=32)
        params["eval_metric"] = "ndcg@10"
        if lambdarank_pair_method is not None:
            params["lambdarank_pair_method"] = lambdarank_pair_method
        if lambdarank_num_pair is not None:
            params["lambdarank_num_pair_per_sample"] = lambdarank_num_pair
    elif objective == "rank:pairwise":
        y = y_raw
    else:
        raise ValueError(objective)
    d = xgb.QuantileDMatrix(X, label=y, group=groups, feature_names=feats)
    return xgb.train(params, d, num_boost_round=rounds)


def _predict(booster, df, feats):
    X = df[feats].to_numpy(np.float32)
    d = xgb.DMatrix(X, feature_names=feats)
    return booster.predict(d)


def _cohort_metrics(preds, y_true, cohort_idx, top_ns=(5, 10)):
    """Per-cohort Spearman, NDCG@K, Precision@K, Top-K mean inhibition.
    Returns BOTH mean and median across cohorts for Spearman.
    """
    out = {"n_cohorts": 0, "Spearman_mean": np.nan, "Spearman_median": np.nan,
           "MAE": np.nan, "RMSE": np.nan}
    for n in top_ns:
        out[f"NDCG@{n}"] = np.nan
        out[f"P@{n}"] = np.nan
        out[f"Top{n}_Inhib"] = np.nan

    mask = np.isfinite(preds) & np.isfinite(y_true)
    if mask.any():
        out["MAE"] = float(mean_absolute_error(y_true[mask], preds[mask]))
        out["RMSE"] = float(np.sqrt(mean_squared_error(y_true[mask], preds[mask])))

    sp = []
    ndcgs = {n: [] for n in top_ns}
    ps = {n: [] for n in top_ns}
    tops = {n: [] for n in top_ns}
    for idxs in cohort_idx:
        t, p = y_true[idxs], preds[idxs]
        finite = np.isfinite(t) & np.isfinite(p)
        if finite.sum() < 2:
            continue
        t, p = t[finite], p[finite]
        if np.ptp(t) == 0 or np.ptp(p) == 0:
            continue
        c, _ = spearmanr(t, p)
        if np.isfinite(c):
            sp.append(c)
        order_pred = np.argsort(-p, kind="stable")
        order_true = np.argsort(-t, kind="stable")
        for n in top_ns:
            if len(t) >= n:
                ndcgs[n].append(ndcg_score(t.reshape(1, -1), p.reshape(1, -1), k=n))
                pred_top = set(order_pred[:n])
                true_top = set(order_true[:n])
                ps[n].append(len(pred_top & true_top) / n)
                tops[n].append(float(np.mean(t[order_pred[:n]])))

    out["n_cohorts"] = len(sp)
    if sp:
        out["Spearman_mean"]   = float(np.mean(sp))
        out["Spearman_median"] = float(np.median(sp))
    for n in top_ns:
        if ndcgs[n]:
            out[f"NDCG@{n}"]      = float(np.mean(ndcgs[n]))
        if ps[n]:
            out[f"P@{n}"]         = float(np.mean(ps[n]))
        if tops[n]:
            out[f"Top{n}_Inhib"]  = float(np.mean(tops[n]))
    return out


def _large_cohort_indices(df, group_cols, min_size=10):
    df = df.reset_index(drop=True)
    return [g.index.values for _, g in df.groupby(group_cols) if len(g) >= min_size]


def _load_data():
    """Read the v3 cache directly + relaxed CSV, compute the 17 new (masked) features."""
    cache_file = ensure_cache("oligo")
    idx = cache_index_col("oligo")
    cache_df = pd.read_parquet(cache_file)
    raw = pd.read_csv(OLIGO_CSV_PROCESSED_AVERAGED)
    common = sorted(set(cache_df.columns) & set(raw.columns) - {idx})
    final_data = pd.merge(cache_df, raw.drop(columns=common), on=idx)

    features_to_ignore = {idx, INHIBITION, "inhibition_percent", "dosage", "split"}
    all_features = sorted(
        set(c for c in final_data.select_dtypes(include=["number"]).columns if c not in features_to_ignore)
        | {VOLUME}
    )
    all_features = [f for f in all_features if f not in COMPETITION]
    final_data = add_new_features(final_data, mask_non_gapmer=True)

    # Use the same NEW feature set the chem-overhaul branch produces in production:
    feats = [f for f in all_features if f not in DROPPED_IN_NEW] + ADDED_IN_NEW
    return final_data, feats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s %(levelname)s | %(message)s",
                        stream=sys.stdout, force=True)

    logger.info("Loading data + features ...")
    final_data, feats = _load_data()
    logger.info("  rows=%d, features=%d (NEW set: cache - {is_hepa,is_all_ps} + 17 masked)",
                len(final_data), len(feats))

    trv = final_data[final_data["split"].isin(["train", "val"])].copy()
    test = final_data[final_data["split"] == "test"].copy()
    y_test = test[INHIBITION].to_numpy(np.float64)
    cust_cohorts = _large_cohort_indices(test, "custom_id", min_size=10)
    cell_cohorts = _large_cohort_indices(test, [CANONICAL_GENE, CELL_LINE], min_size=10)
    logger.info("  train+val=%d, test=%d  |  custom_id cohorts(>=10)=%d  |  gene x cell cohorts(>=10)=%d",
                len(trv), len(test), len(cust_cohorts), len(cell_cohorts))

    # Variant grid: (label, objective, group_choice, rounds, pair_method, num_pair)
    variants = [
        ("L2",                       "reg:squarederror", None,        700,  None,   None),
        ("pair_custid",              "rank:pairwise",    "custom_id", 700,  None,   None),
        ("pair_gene_cell",           "rank:pairwise",    "gene_cell", 700,  None,   None),
        ("pair_custid_long",         "rank:pairwise",    "custom_id", 1400, None,   None),
        ("ndcg_custid_default",      "rank:ndcg",        "custom_id", 700,  None,   None),
        ("ndcg_custid_mean_8",       "rank:ndcg",        "custom_id", 700,  "mean", 8),
        ("ndcg_custid_topk_4",       "rank:ndcg",        "custom_id", 700,  "topk", 4),
        ("ndcg_gene_cell_default",   "rank:ndcg",        "gene_cell", 700,  None,   None),
        ("ndcg_gene_cell_mean_8",    "rank:ndcg",        "gene_cell", 700,  "mean", 8),
    ]

    rows = []
    total = len(variants) * len(SEEDS)
    done = 0
    t0 = time.perf_counter()

    for label, obj, gc, rounds, pair_method, num_pair in variants:
        for seed in SEEDS:
            done += 1
            booster = _train(
                trv, feats, seed=seed, objective=obj, group_choice=gc, rounds=rounds,
                lambdarank_pair_method=pair_method, lambdarank_num_pair=num_pair,
            )
            preds = _predict(booster, test, feats)
            m_cust = _cohort_metrics(preds, y_test, cust_cohorts)
            m_cell = _cohort_metrics(preds, y_test, cell_cohorts)
            elapsed = time.perf_counter() - t0
            eta = elapsed / done * (total - done)
            logger.info("[%2d/%d eta=%.0fs] %s seed=%d  "
                        "cust_sp_mean=%.4f  cust_sp_med=%.4f  cell_sp_mean=%.4f  cust_ndcg10=%.4f",
                        done, total, eta, label, seed,
                        m_cust["Spearman_mean"], m_cust["Spearman_median"],
                        m_cell["Spearman_mean"], m_cust["NDCG@10"])
            rows.append(dict(
                variant=label, objective=obj, group_choice=gc or "-", rounds=rounds,
                pair_method=pair_method or "-", num_pair=num_pair if num_pair is not None else 0,
                seed=seed,
                **{f"cust_{k}": v for k, v in m_cust.items()},
                **{f"cell_{k}": v for k, v in m_cell.items()},
            ))

    long = pd.DataFrame(rows)
    out_csv = Path(__file__).parent / "_lambdarank_relaxed_results.csv"
    long.to_csv(out_csv, index=False)
    logger.info("Saved long results to %s", out_csv)

    # Aggregate across seeds. For each cohort grouping, report mean+std of the seed-mean.
    agg_cols = {
        "cust_Spearman_mean":   ("cust_Spearman_mean", "mean"),
        "cust_Spearman_med":    ("cust_Spearman_median", "mean"),
        "cust_NDCG_5":          ("cust_NDCG@5", "mean"),
        "cust_NDCG_10":         ("cust_NDCG@10", "mean"),
        "cust_P_5":             ("cust_P@5", "mean"),
        "cust_P_10":            ("cust_P@10", "mean"),
        "cust_Top5":            ("cust_Top5_Inhib", "mean"),
        "cust_Top10":           ("cust_Top10_Inhib", "mean"),
        "cell_Spearman_mean":   ("cell_Spearman_mean", "mean"),
        "cell_Spearman_med":    ("cell_Spearman_median", "mean"),
        "cell_NDCG_5":          ("cell_NDCG@5", "mean"),
        "cell_NDCG_10":         ("cell_NDCG@10", "mean"),
        "cell_P_5":             ("cell_P@5", "mean"),
        "cell_P_10":            ("cell_P@10", "mean"),
        "cell_Top5":            ("cell_Top5_Inhib", "mean"),
        "cell_Top10":           ("cell_Top10_Inhib", "mean"),
        "MAE":                  ("cust_MAE", "mean"),
    }
    grid = long.groupby("variant", as_index=False).agg(**agg_cols)

    pd.set_option("display.width", 240)
    pd.set_option("display.float_format", lambda x: f"{x:+.4f}" if pd.notnull(x) else "    nan")

    print(f"\n=== Custom-id cohort metrics ({len(SEEDS)} seeds, mean across seeds) ===")
    print(grid[["variant",
                "cust_Spearman_mean", "cust_Spearman_med",
                "cust_NDCG_5", "cust_NDCG_10", "cust_P_5", "cust_P_10",
                "cust_Top5", "cust_Top10"]].to_string(index=False))

    print(f"\n=== Gene x cell-line cohort metrics ({len(SEEDS)} seeds, mean across seeds) ===")
    print(grid[["variant",
                "cell_Spearman_mean", "cell_Spearman_med",
                "cell_NDCG_5", "cell_NDCG_10", "cell_P_5", "cell_P_10",
                "cell_Top5", "cell_Top10"]].to_string(index=False))

    print(f"\n=== Magnitude (MAE only meaningful for L2) ===")
    print(grid[["variant", "MAE"]].to_string(index=False))


if __name__ == "__main__":
    main()
