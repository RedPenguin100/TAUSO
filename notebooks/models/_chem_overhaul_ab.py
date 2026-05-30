"""A/B between the old chem feature set and the new one shipped in feat/expr-prefix-and-chem-review.

Reads ONLY the v3 Zenodo wide cache (oligo_full_combined_v3.parquet) -- the same source
that train_model.py + evaluate_all.py treat as canonical. Deliberately ignores the loose
parquet shards in .tauso_data/features/oligo/, which are leftovers from past pipeline
runs and would otherwise double-count features.

OLD set = the v3 cache columns as-is, including the now-dropped `is_hepa` and `is_all_ps`.
NEW set = OLD minus {is_hepa, is_all_ps} plus the 17 new chem/architecture features
added across this branch. The 17 are computed in-memory from CHEMICAL_PATTERN /
PS_PATTERN / SEQUENCE / MODIFICATION on the source CSV.

The expression-feature rename (`target_expression` -> `expr_target` etc.) is information-
neutral for an XGBoost model trained from a numpy matrix; the rename is not part of the A/B.

Trains 3 seeds per variant with SANE hyperparams and reports cust/gene-cell Spearman,
NDCG@10, Top10_Inhib on the test split. COMPETITION outputs are excluded from features.
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
from notebooks.models.utility import _apply_other_as_nan
from tauso.common.modifications import get_longest_dna_gap
from tauso.data.consts import CANONICAL_GENE, CELL_LINE, CHEMICAL_PATTERN, INHIBITION, PS_PATTERN, SEQUENCE, VOLUME
from tauso.features.sequence.toxicity_features import cpg_ps_count
from tauso.populate.feature_cache import cache_path, ensure_cache, index_col as cache_index_col

logger = logging.getLogger(__name__)

SANE_BASE = {
    "tree_method": "hist", "device": "cuda",
    "max_depth": 6, "learning_rate": 0.05,
    "subsample": 0.8, "colsample_bytree": 0.6, "min_child_weight": 20,
    "reg_lambda": 5.0, "reg_alpha": 1.0,
    "objective": "reg:squarederror",
}
ROUNDS = 700
SEEDS = (1, 2, 3)

DROPPED_IN_NEW = ["is_hepa", "is_all_ps"]
ADDED_IN_NEW = [
    "chem_1st_gen",
    "arch_wing5_len", "arch_gap_len", "arch_wing3_len", "arch_is_gapmer",
    "hybr_dna_rna_wing_dg_imbalance",
    "hybr_dna_rna_wing5_minus_gap_dg",
    "hybr_dna_rna_wing3_minus_gap_dg",
    "ps_wing5_count", "ps_gap_count", "ps_wing3_count",
    "frac_mod_with_ps", "frac_dna_with_ps",
    "seq5_is_purine", "seq5_is_g", "seq5_is_t",
    "tox_cpg_ps_count",
]


def _arch_lens(pat):
    if not isinstance(pat, str) or not pat:
        return 0, 0, 0
    start, end, gap_len = get_longest_dna_gap(pat)
    if gap_len == 0:
        return 0, 0, 0
    return start, gap_len, len(pat) - end


def _ps_placement(chem_pat, ps_pat):
    if not isinstance(chem_pat, str) or not isinstance(ps_pat, str):
        return 0, 0, 0
    gap_start, gap_end, gap_len = get_longest_dna_gap(chem_pat)
    if gap_len == 0:
        return 0, 0, 0
    return (
        ps_pat[:gap_start].count("*"),
        ps_pat[gap_start:gap_end].count("*"),
        ps_pat[gap_end:].count("*"),
    )


def _frac_with_ps(chem_pat, ps_pat, marker_chars):
    if not isinstance(chem_pat, str) or not isinstance(ps_pat, str):
        return 0.0
    n = min(len(chem_pat), len(ps_pat))
    qualifying, ps = 0, 0
    for i in range(n):
        if chem_pat[i] in marker_chars:
            qualifying += 1
            if ps_pat[i] == "*":
                ps += 1
    return ps / qualifying if qualifying else 0.0


def add_new_features(df):
    """Adds the 17 new features in-place. Idempotent."""
    from tauso.data.consts import MODIFICATION

    out = df.copy()

    # 1. chem_1st_gen
    has_moe = out[MODIFICATION].str.contains("MOE", na=False)
    has_high_aff = out[MODIFICATION].str.contains("LNA|cEt", na=False)
    out["chem_1st_gen"] = (~(has_moe | has_high_aff)).astype(int)

    # 2. arch_* — gapmer geometry
    arch_tuples = out[CHEMICAL_PATTERN].map(_arch_lens)
    out["arch_wing5_len"] = arch_tuples.map(lambda t: t[0]).astype(int)
    out["arch_gap_len"]   = arch_tuples.map(lambda t: t[1]).astype(int)
    out["arch_wing3_len"] = arch_tuples.map(lambda t: t[2]).astype(int)
    out["arch_is_gapmer"] = arch_tuples.map(
        lambda t: 1 if t[0] > 0 and t[1] > 0 and t[2] > 0 else 0
    ).astype(int)

    # 3. hybr wing-gap asymmetry (derived from existing columns)
    out["hybr_dna_rna_wing_dg_imbalance"] = out["hybr_dna_rna_dg_wing5"] - out["hybr_dna_rna_dg_wing3"]
    out["hybr_dna_rna_wing5_minus_gap_dg"] = out["hybr_dna_rna_dg_wing5"] - out["hybr_dna_rna_dg_gap"]
    out["hybr_dna_rna_wing3_minus_gap_dg"] = out["hybr_dna_rna_dg_wing3"] - out["hybr_dna_rna_dg_gap"]

    # 4. PS placement
    triples = out.apply(lambda r: _ps_placement(r[CHEMICAL_PATTERN], r[PS_PATTERN]), axis=1)
    out["ps_wing5_count"] = triples.map(lambda t: t[0]).astype(int)
    out["ps_gap_count"]   = triples.map(lambda t: t[1]).astype(int)
    out["ps_wing3_count"] = triples.map(lambda t: t[2]).astype(int)

    # 5. mod x backbone interaction
    out["frac_mod_with_ps"] = out.apply(
        lambda r: _frac_with_ps(r[CHEMICAL_PATTERN], r[PS_PATTERN], {"M", "C", "L"}), axis=1
    ).astype(float)
    out["frac_dna_with_ps"] = out.apply(
        lambda r: _frac_with_ps(r[CHEMICAL_PATTERN], r[PS_PATTERN], {"d"}), axis=1
    ).astype(float)

    # 6. 5'-terminal base
    seq_5p = out[SEQUENCE].str[0]
    out["seq5_is_purine"] = seq_5p.isin(["A", "G"]).astype(int)
    out["seq5_is_g"]      = (seq_5p == "G").astype(int)
    out["seq5_is_t"]      = seq_5p.isin(["T", "U"]).astype(int)

    # 7. tox_cpg_ps_count
    out["tox_cpg_ps_count"] = out.apply(
        lambda r: cpg_ps_count(r[SEQUENCE], r[PS_PATTERN]), axis=1
    ).astype(int)

    return out


def _train(trv, feats, seed):
    y = trv[INHIBITION].to_numpy(np.float64)
    X = trv[feats].to_numpy(np.float32)
    params = {**SANE_BASE, "seed": seed}
    d = xgb.QuantileDMatrix(X, label=y, feature_names=feats)
    return xgb.train(params, d, num_boost_round=ROUNDS)


def _predict(booster, df, feats):
    X = df[feats].to_numpy(np.float32)
    d = xgb.DMatrix(X, feature_names=feats)
    return booster.predict(d)


def _cohort_metrics(preds, y_true, cohort_idx, top_n=(10,)):
    out = {"n_cohorts": 0, "Spearman": np.nan, "MAE": np.nan, "RMSE": np.nan}
    for n in top_n:
        out[f"NDCG@{n}"] = np.nan
        out[f"Top{n}_Inhib"] = np.nan
    mask = np.isfinite(preds) & np.isfinite(y_true)
    if mask.any():
        out["MAE"] = float(mean_absolute_error(y_true[mask], preds[mask]))
        out["RMSE"] = float(np.sqrt(mean_squared_error(y_true[mask], preds[mask])))
    sp, ndcgs, tops = [], {n: [] for n in top_n}, {n: [] for n in top_n}
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
        order = np.argsort(-p, kind="stable")
        for n in top_n:
            if len(t) >= n:
                ndcgs[n].append(ndcg_score(t.reshape(1, -1), p.reshape(1, -1), k=n))
                tops[n].append(float(np.mean(t[order[:n]])))
    out["n_cohorts"] = len(sp)
    if sp:
        out["Spearman"] = float(np.mean(sp))
    for n in top_n:
        if ndcgs[n]:
            out[f"NDCG@{n}"] = float(np.mean(ndcgs[n]))
        if tops[n]:
            out[f"Top{n}_Inhib"] = float(np.mean(tops[n]))
    return out


def _large_cohort_indices(df, group_cols, min_size=10):
    df = df.reset_index(drop=True)
    return [g.index.values for _, g in df.groupby(group_cols) if len(g) >= min_size]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s %(levelname)s | %(message)s",
                        stream=sys.stdout, force=True)

    logger.info("Loading the v3 wide cache ONLY (ignoring loose shards) ...")
    cache_file = ensure_cache("oligo")
    idx = cache_index_col("oligo")
    cache_features = pd.read_parquet(cache_file)
    logger.info("  cache: %s -> %d rows, %d cols", cache_file, len(cache_features), cache_features.shape[1])

    logger.info("Loading the source CSV (raw metadata: sequence, patterns, split, cohorts) ...")
    raw = pd.read_csv(OLIGO_CSV_PROCESSED_AVERAGED)
    logger.info("  source: %d rows, %d cols", len(raw), raw.shape[1])

    # Same merge-dedup pattern as load_and_validate_final_data: drop columns common to both
    # the cache and the source from the source side, then merge on the index.
    common = sorted(set(cache_features.columns) & set(raw.columns) - {idx})
    logger.info("  dropping %d shared cols from source: %s", len(common), common[:10])
    final_data = pd.merge(cache_features, raw.drop(columns=common), on=idx)
    logger.info("  merged final_data: %d rows, %d cols", len(final_data), final_data.shape[1])

    # Trainable feature list = numeric columns minus the things that aren't features.
    features_to_ignore = {idx, INHIBITION, "inhibition_percent", "dosage", "split"}
    all_features = sorted(
        set(c for c in final_data.select_dtypes(include=["number"]).columns if c not in features_to_ignore)
        | {VOLUME}
    )
    # COMPETITION outputs are baselines, not inputs.
    all_features = [f for f in all_features if f not in COMPETITION]
    logger.info("  trainable features (cache - COMPETITION): %d", len(all_features))

    # Compute the 17 new features in-memory from the source columns. These columns are
    # already on final_data (CHEMICAL_PATTERN, PS_PATTERN, SEQUENCE, MODIFICATION).
    logger.info("Adding the 17 new chem features in-memory ...")
    t0 = time.perf_counter()
    final_data = add_new_features(final_data)
    logger.info("  done in %.1fs", time.perf_counter() - t0)

    # Sanity check
    for f in ADDED_IN_NEW:
        assert f in final_data.columns, f"missing new feature: {f}"

    # OLD feature set: the cache feature list as-is (already includes is_hepa/is_all_ps).
    feats_old = list(all_features)

    # NEW feature set: drop the two retired flags, add the 17 new ones (after deduping).
    feats_new = [f for f in all_features if f not in DROPPED_IN_NEW]
    # The new features are NOT in `all_features` from load_and_validate_final_data (they
    # were just added by add_new_features), so it's safe to extend.
    feats_new = feats_new + ADDED_IN_NEW

    logger.info("  OLD set: %d features  | NEW set: %d features (diff: %+d)",
                len(feats_old), len(feats_new), len(feats_new) - len(feats_old))

    trv = final_data[final_data["split"].isin(["train", "val"])].copy()
    test = final_data[final_data["split"] == "test"].copy()
    y_test = test[INHIBITION].to_numpy(np.float64)
    cust_cohorts = _large_cohort_indices(test, "custom_id", min_size=10)
    cell_cohorts = _large_cohort_indices(test, [CANONICAL_GENE, CELL_LINE], min_size=10)
    logger.info("  train+val=%d, test=%d, custom_id cohorts(>=10)=%d, gene x cell cohorts(>=10)=%d",
                len(trv), len(test), len(cust_cohorts), len(cell_cohorts))

    variants = (("OLD", feats_old), ("NEW", feats_new))
    rows = []
    total = len(variants) * len(SEEDS)
    done = 0
    t0 = time.perf_counter()

    for label, feats in variants:
        for seed in SEEDS:
            done += 1
            booster = _train(trv, feats, seed)
            preds = _predict(booster, test, feats)
            m_cust = _cohort_metrics(preds, y_test, cust_cohorts)
            m_cell = _cohort_metrics(preds, y_test, cell_cohorts)
            elapsed = time.perf_counter() - t0
            eta = elapsed / done * (total - done)
            logger.info("[%2d/%d eta=%.0fs] variant=%s seed=%d nfeats=%d  "
                        "cust_sp=%.4f  cell_sp=%.4f  cust_ndcg10=%.4f",
                        done, total, eta, label, seed, len(feats),
                        m_cust["Spearman"], m_cell["Spearman"], m_cust["NDCG@10"])
            rows.append(dict(variant=label, seed=seed, n_features=len(feats),
                             **{f"cust_{k}": v for k, v in m_cust.items()},
                             **{f"cell_{k}": v for k, v in m_cell.items()}))

    long = pd.DataFrame(rows)
    out_csv = Path(__file__).parent / "_chem_overhaul_ab_results.csv"
    long.to_csv(out_csv, index=False)
    logger.info("Saved long results to %s", out_csv)

    grid = long.groupby(["variant"], as_index=False).agg(
        sp_custom_id_mean   =("cust_Spearman", "mean"),
        sp_custom_id_std    =("cust_Spearman", "std"),
        sp_gene_x_cell_mean =("cell_Spearman", "mean"),
        sp_gene_x_cell_std  =("cell_Spearman", "std"),
        ndcg10_custom_id    =("cust_NDCG@10", "mean"),
        ndcg10_gene_x_cell  =("cell_NDCG@10", "mean"),
        top10_inhib_custom_id   =("cust_Top10_Inhib", "mean"),
        top10_inhib_gene_x_cell =("cell_Top10_Inhib", "mean"),
        mae_custom_id       =("cust_MAE", "mean"),
        n_features          =("n_features", "max"),
    )
    pd.set_option("display.width", 220)
    pd.set_option("display.float_format", lambda x: f"{x:+.4f}" if pd.notnull(x) else "    nan")
    print(f"\n=== Old vs New chem feature set ({len(SEEDS)} seeds each, SANE XGBoost, test split) ===")
    print(grid.to_string(index=False))

    # OldNew delta for the eye
    if {"OLD", "NEW"}.issubset(set(long["variant"])):
        a = long[long["variant"] == "OLD"].mean(numeric_only=True)
        b = long[long["variant"] == "NEW"].mean(numeric_only=True)
        print("\nNEW − OLD (positive = new set wins):")
        for k in ("cust_Spearman", "cell_Spearman", "cust_NDCG@10", "cell_NDCG@10",
                  "cust_Top10_Inhib", "cell_Top10_Inhib", "cust_MAE"):
            print(f"  {k:20s}  {b[k] - a[k]:+.4f}")

    if "oligo_ai_score" in test.columns:
        oa = test["oligo_ai_score"].to_numpy(np.float64)
        m_cust = _cohort_metrics(oa, y_test, cust_cohorts)
        m_cell = _cohort_metrics(oa, y_test, cell_cohorts)
        print(f"\n=== OligoAI competitor reference (no training) ===")
        print(f"  Spearman/custom_id   = {m_cust['Spearman']:+.4f}   NDCG@10 = {m_cust['NDCG@10']:+.4f}   Top10 = {m_cust['Top10_Inhib']:+.2f}")
        print(f"  Spearman/gene×cell   = {m_cell['Spearman']:+.4f}   NDCG@10 = {m_cell['NDCG@10']:+.4f}   Top10 = {m_cell['Top10_Inhib']:+.2f}")


if __name__ == "__main__":
    main()
