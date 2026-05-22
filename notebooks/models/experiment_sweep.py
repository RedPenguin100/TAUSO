"""
Overnight experiment sweep: compare model configurations for ASO efficacy prediction.

Varies: loss function, feature importance method, cohort weighting, evaluation
grouping strategy, data size, and feature selection approach.

Each result is saved incrementally after every experiment — partial results
survive a crash or early termination.

Usage:
    python -u -m notebooks.models.experiment_sweep 2>&1 | tee sweep.log
    python -u -m notebooks.models.experiment_sweep --cpus 8 --output results/sweep.json

Runtime: ~3–6 hours on 8 CPUs. Memory: <4 GB peak.
"""
import argparse
import gc
import json
import logging
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional

import numpy as np
import optuna
import pandas as pd
import xgboost as xgb
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error

from notebooks.models.utility import load_and_validate_final_data
from tauso.data.consts import CANONICAL_GENE, CELL_LINE, INHIBITION

logger = logging.getLogger(__name__)
optuna.logging.set_verbosity(optuna.logging.WARNING)

# ---------------------------------------------------------------------------
# Experiment configs
# ---------------------------------------------------------------------------

@dataclass
class ExperimentConfig:
    name: str
    n_train: int          # rows sampled from train split
    loss: str             # "L2", "huber"
    importance: str       # "gain", "permutation"
    weighting: str        # "uniform", "inv_cohort"
    eval_cols: List[str]  # columns to group by for eval
    selection: str        # "backward", "threshold"
    n_optuna: int = 20


OBJECTIVES = {
    "L2":    "reg:squarederror",
    "huber": "reg:pseudohubererror",
}

GC = [CANONICAL_GENE, CELL_LINE]   # gene × cell line
GO = [CANONICAL_GENE]              # gene only
CO = [CELL_LINE]                   # cell line only
CI = ["custom_id"]                 # per-experiment batch

EXPERIMENTS = [
    # ---------- baseline: vary only data size ----------
    ExperimentConfig("baseline_4k",    4000, "L2",    "gain",        "uniform",    GC, "backward"),
    ExperimentConfig("baseline_6k",    6000, "L2",    "gain",        "uniform",    GC, "backward"),
    # ---------- loss function ----------
    ExperimentConfig("huber_4k",       4000, "huber", "gain",        "uniform",    GC, "backward"),
    ExperimentConfig("huber_6k",       6000, "huber", "gain",        "uniform",    GC, "backward"),
    # ---------- feature importance method ----------
    ExperimentConfig("perm_4k",        4000, "L2",    "permutation", "uniform",    GC, "backward"),
    ExperimentConfig("perm_6k",        6000, "L2",    "permutation", "uniform",    GC, "backward"),
    # ---------- cohort weighting ----------
    ExperimentConfig("inv_wt_4k",      4000, "L2",    "gain",        "inv_cohort", GC, "backward"),
    ExperimentConfig("inv_wt_6k",      6000, "L2",    "gain",        "inv_cohort", GC, "backward"),
    # ---------- evaluation grouping ----------
    ExperimentConfig("gene_only_4k",   4000, "L2",    "gain",        "uniform",    GO, "backward"),
    ExperimentConfig("cell_only_4k",   4000, "L2",    "gain",        "uniform",    CO, "backward"),
    ExperimentConfig("custom_id_4k",   4000, "L2",    "gain",        "uniform",    CI, "backward"),
    # ---------- feature selection strategy ----------
    ExperimentConfig("threshold_6k",   6000, "L2",    "gain",        "uniform",    GC, "threshold"),
    # ---------- combinations ----------
    ExperimentConfig("huber_inv_6k",   6000, "huber", "gain",        "inv_cohort", GC, "backward"),
    ExperimentConfig("perm_inv_4k",    4000, "L2",    "permutation", "inv_cohort", GC, "backward"),
    ExperimentConfig("huber_perm_4k",  4000, "huber", "permutation", "uniform",    GC, "backward"),
]

# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

def _cohort_indices(df, eval_cols, min_size=20):
    df = df.reset_index(drop=True)
    return [
        g.index.values
        for _, g in df.groupby(eval_cols)
        if len(g) >= min_size
    ]


def _cohort_spearman(y_true, preds, indices):
    scores = [
        spearmanr(y_true[idx], preds[idx])[0]
        for idx in indices
        if len(idx) >= 3
    ]
    scores = [s for s in scores if not np.isnan(s)]
    return float(np.median(scores)) if scores else float("nan"), len(scores)


def _cohort_weights(df, eval_cols):
    sizes = df.groupby(eval_cols)[eval_cols[0]].transform("count")
    w = 1.0 / sizes.values.astype(np.float32)
    return w / w.mean()


def _sample_train(final_data, n_rows, seed=1):
    """Sample n_rows from train split, capped per cohort for diversity."""
    train = final_data[final_data["split"] == "train"]
    n_cohorts = train.groupby(GC).ngroups
    per_cohort = max(5, n_rows // n_cohorts)
    rows = []
    for _, g in train.groupby(GC):
        rows.append(g.sample(min(per_cohort, len(g)), random_state=seed))
    sampled = pd.concat(rows)
    if len(sampled) > n_rows:
        sampled = sampled.sample(n_rows, random_state=seed)
    return sampled.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Optuna tuning (light, per experiment)
# ---------------------------------------------------------------------------

def _tune(X_tr, y_tr, X_vl, y_vl, val_optuna_groups, objective, n_trials, nthread, seed, weights_tr=None):
    dtrain = xgb.DMatrix(X_tr, label=y_tr)
    if weights_tr is not None:
        dtrain.set_weight(weights_tr)
    dval = xgb.DMatrix(X_vl, label=y_vl)

    def fast_spearman(preds):
        # idx = global indices into val; y = already-sliced values for that group
        scores = [
            spearmanr(y, preds[idx])[0]
            for idx, y in val_optuna_groups
            if not np.isnan(spearmanr(y, preds[idx])[0])
        ]
        return float(np.nanmean(scores)) if scores else 0.0

    def objective_fn(trial):
        params = {
            "tree_method": "hist", "device": "cpu", "nthread": nthread,
            "objective": objective, "seed": seed,
            "max_depth":         trial.suggest_int("max_depth", 2, 6),
            "learning_rate":     trial.suggest_float("learning_rate", 0.01, 0.15, log=True),
            "subsample":         trial.suggest_float("subsample", 0.4, 0.9),
            "colsample_bytree":  trial.suggest_float("colsample_bytree", 0.3, 0.9),
            "min_child_weight":  trial.suggest_int("min_child_weight", 10, 200),
            "gamma":             trial.suggest_float("gamma", 0.01, 5.0, log=True),
            "reg_alpha":         trial.suggest_float("reg_alpha", 0.01, 50.0, log=True),
            "reg_lambda":        trial.suggest_float("reg_lambda", 0.01, 50.0, log=True),
        }
        if objective == "reg:pseudohubererror":
            params["huber_slope"] = trial.suggest_float("huber_slope", 0.5, 10.0, log=True)

        bst = xgb.train(
            params, dtrain,
            num_boost_round=trial.suggest_int("num_boost_round", 200, 1000, step=100),
            evals=[(dval, "val")], early_stopping_rounds=30, verbose_eval=False,
        )
        trial.set_user_attr("best_iteration", bst.best_iteration)
        return fast_spearman(bst.predict(dval))

    study = optuna.create_study(direction="maximize",
                                sampler=optuna.samplers.TPESampler(seed=seed))
    study.optimize(objective_fn, n_trials=n_trials, show_progress_bar=False)

    best = study.best_params.copy()
    best["num_boost_round"] = study.best_trial.user_attrs["best_iteration"]
    best.update({"tree_method": "hist", "device": "cpu", "nthread": nthread,
                 "objective": objective, "seed": seed})
    logger.info("  Optuna done | best_spearman=%.4f | rounds=%d depth=%d lr=%.4f",
                study.best_value, best["num_boost_round"],
                best["max_depth"], best["learning_rate"])
    return best


# ---------------------------------------------------------------------------
# Feature importance
# ---------------------------------------------------------------------------

def _gain_importance(bst, features):
    raw = bst.get_score(importance_type="gain")
    return {f: raw.get(f, 0.0) for f in features}


def _permutation_importance(bst, X_val, y_val, val_select_idx, features, seed=0):
    rng = np.random.default_rng(seed)
    base_preds = bst.inplace_predict(X_val)
    if hasattr(base_preds, "get"):
        base_preds = base_preds.get()
    baseline, _ = _cohort_spearman(y_val, base_preds, val_select_idx)
    importances = {}
    for i, feat in enumerate(features):
        X_p = X_val.copy()
        X_p[:, i] = rng.permutation(X_p[:, i])
        preds = bst.inplace_predict(X_p)
        if hasattr(preds, "get"):
            preds = preds.get()
        perturbed, _ = _cohort_spearman(y_val, preds, val_select_idx)
        importances[feat] = baseline - perturbed  # higher = more important
    return importances


# ---------------------------------------------------------------------------
# Backward RFE
# ---------------------------------------------------------------------------

PARSIMONY = 0.005
EARLY_STOP_DROP = 0.15


def _backward_rfe(X_tr, y_tr, X_vl, y_vl, features, params, val_eval_idx,
                  val_select_idx, importance_type, weights_tr, seed):
    feat_to_idx = {f: i for i, f in enumerate(features)}
    current = list(features)
    history = []
    peak = -1.0
    num_rounds = params.pop("num_boost_round", 500)
    p = params.copy()

    logger.info("  RFE start | features=%d importance=%s", len(current), importance_type)

    while len(current) > 1:
        idxs = [feat_to_idx[f] for f in current]
        Xtr_c, Xvl_c = X_tr[:, idxs], X_vl[:, idxs]

        dtr = xgb.QuantileDMatrix(Xtr_c, label=y_tr, feature_names=current)
        if weights_tr is not None:
            dtr = xgb.DMatrix(Xtr_c, label=y_tr, feature_names=current)
            dtr.set_weight(weights_tr)

        bst = xgb.train(p, dtr, num_boost_round=num_rounds, verbose_eval=False)

        preds = bst.inplace_predict(Xvl_c)
        if hasattr(preds, "get"):
            preds = preds.get()

        val_spear, _ = _cohort_spearman(y_vl, preds, val_select_idx)
        peak = max(peak, val_spear)

        history.append({
            "n_features": len(current),
            "val_spearman": val_spear,
        })

        if val_spear < peak - EARLY_STOP_DROP:
            logger.info("  RFE early stop | n=%d val=%.4f peak=%.4f", len(current), val_spear, peak)
            break

        if importance_type == "permutation":
            imp = _permutation_importance(bst, Xvl_c, y_vl, val_select_idx, current, seed=seed)
        else:
            imp = _gain_importance(bst, current)

        ranked = sorted(imp.items(), key=lambda x: x[1])
        drop_n = min(5 if len(current) > 100 else 1, len(current) - 1)
        to_drop = {f for f, _ in ranked[:drop_n]}
        current = [f for f in current if f not in to_drop]

        del dtr, bst, Xtr_c, Xvl_c
        gc.collect()

    # parsimony: pick smallest set within PARSIMONY of peak
    df_hist = pd.DataFrame(history)
    threshold = df_hist["val_spearman"].max() - PARSIMONY
    best_row = df_hist[df_hist["val_spearman"] >= threshold].iloc[-1]
    n_selected = int(best_row["n_features"])
    best_val = float(best_row["val_spearman"])

    # find the feature list at that step
    total_steps = len(df_hist)
    step_idx = df_hist[df_hist["val_spearman"] >= threshold].index[-1]
    # reconstruct which features were current at that step
    # (we track n_features but not the list — report n only)
    logger.info("  RFE done | selected=%d best_val_spearman=%.4f", n_selected, best_val)
    return n_selected, best_val, df_hist


# ---------------------------------------------------------------------------
# Threshold sweep (alternative to RFE)
# ---------------------------------------------------------------------------

THRESHOLDS = [10, 25, 50, 75, 100, 150, 200, 300, 500, 700, None]  # None = all


def _threshold_sweep(X_tr, y_tr, X_vl, y_vl, features, params, val_eval_idx, weights_tr, seed):
    num_rounds = params.pop("num_boost_round", 500)
    p = params.copy()

    # Train full model to get importance ranking
    dtr_full = xgb.DMatrix(X_tr, label=y_tr, feature_names=features)
    if weights_tr is not None:
        dtr_full.set_weight(weights_tr)
    bst_full = xgb.train(p, dtr_full, num_boost_round=num_rounds, verbose_eval=False)
    raw_imp = bst_full.get_score(importance_type="gain")
    ranked_feats = sorted(features, key=lambda f: raw_imp.get(f, 0.0), reverse=True)

    results = []
    for k in THRESHOLDS:
        top_k = ranked_feats[:k] if k is not None else ranked_feats
        idxs = [features.index(f) for f in top_k]
        Xtr_k, Xvl_k = X_tr[:, idxs], X_vl[:, idxs]

        dtr_k = xgb.DMatrix(Xtr_k, label=y_tr, feature_names=top_k)
        if weights_tr is not None:
            dtr_k.set_weight(weights_tr)
        bst_k = xgb.train(p, dtr_k, num_boost_round=num_rounds, verbose_eval=False)

        preds = bst_k.predict(xgb.DMatrix(Xvl_k, feature_names=top_k))
        spear, n_cohorts = _cohort_spearman(y_vl, preds, val_eval_idx)
        mae = float(mean_absolute_error(y_vl, preds))
        results.append({"k": k or len(features), "val_spearman": spear, "mae": mae})
        logger.info("  Threshold k=%-4s | spearman=%.4f | mae=%.3f | cohorts=%d",
                    str(k) if k else "all", spear, mae, n_cohorts)

    best = max(results, key=lambda r: r["val_spearman"])
    return best["k"], best["val_spearman"], results


# ---------------------------------------------------------------------------
# Single experiment runner
# ---------------------------------------------------------------------------

def run_experiment(cfg: ExperimentConfig, final_data, features, all_val_df, nthread, seed=1):
    logger.info("\n" + "="*60)
    logger.info("EXPERIMENT: %s", cfg.name)
    logger.info("  loss=%s importance=%s weighting=%s eval=%s selection=%s n_train=%d",
                cfg.loss, cfg.importance, cfg.weighting, cfg.eval_cols, cfg.selection, cfg.n_train)

    t0 = time.time()

    # --- sample train ---
    train_df = _sample_train(final_data, cfg.n_train, seed=seed)
    val_df = all_val_df.copy()
    logger.info("  train=%d val=%d", len(train_df), len(val_df))

    X_tr = train_df[features].values.astype(np.float32)
    y_tr = train_df[INHIBITION].values.astype(np.float32)
    X_vl = val_df[features].values.astype(np.float32)
    y_vl = val_df[INHIBITION].values.astype(np.float32)

    # --- cohort indices for evaluation ---
    min_size = 10 if cfg.eval_cols == ["custom_id"] else 20
    val_eval_idx = _cohort_indices(val_df, cfg.eval_cols, min_size=min_size)
    logger.info("  val eval cohorts (>=%d rows): %d", min_size, len(val_eval_idx))

    # --- sample weights ---
    weights_tr = _cohort_weights(train_df, cfg.eval_cols) if cfg.weighting == "inv_cohort" else None

    # --- Optuna val groups ---
    val_optuna_groups = [
        (g.index.values, y_vl[g.index.values])
        for _, g in val_df.reset_index(drop=True).groupby(cfg.eval_cols)
        if len(g) >= 20
    ]
    logger.info("  Optuna val groups: %d | running %d trials...", len(val_optuna_groups), cfg.n_optuna)

    objective = OBJECTIVES[cfg.loss]

    best_params = _tune(
        X_tr, y_tr, X_vl, y_vl,
        val_optuna_groups, objective,
        n_trials=cfg.n_optuna, nthread=nthread, seed=seed,
        weights_tr=weights_tr,
    )

    # --- feature selection ---
    if cfg.selection == "backward":
        n_selected, best_spear, history = _backward_rfe(
            X_tr, y_tr, X_vl, y_vl, features, best_params, val_eval_idx,
            val_eval_idx, cfg.importance, weights_tr, seed=seed,
        )
        extra = {"rfe_history_len": len(history)}
    else:
        n_selected, best_spear, sweep_results = _threshold_sweep(
            X_tr, y_tr, X_vl, y_vl, features, best_params, val_eval_idx,
            weights_tr, seed=seed,
        )
        extra = {"threshold_curve": sweep_results}

    elapsed = time.time() - t0

    result = {
        "name":         cfg.name,
        "config":       asdict(cfg),
        "val_spearman": best_spear,
        "n_features":   n_selected,
        "elapsed_s":    round(elapsed, 1),
        **extra,
    }

    logger.info("  RESULT | spearman=%.4f | features=%d | time=%.0fs",
                best_spear, n_selected, elapsed)

    del X_tr, X_vl, y_tr, y_vl, train_df
    gc.collect()

    return result


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpus",   type=int, default=8)
    parser.add_argument("--seed",   type=int, default=1)
    parser.add_argument("--output", type=Path,
                        default=Path("notebooks/models/SeenOligoModel/sweep_results.json"))
    parser.add_argument("--resume", action="store_true",
                        help="Skip experiments whose name is already in the output file.")
    parser.add_argument("--only",   nargs="*", default=None,
                        help="Run only these experiment names.")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING"])
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s | %(message)s",
        stream=sys.stdout, force=True,
    )

    # Load results so far (for --resume)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    completed = {}
    if args.resume and args.output.exists():
        with open(args.output) as f:
            for r in json.load(f):
                completed[r["name"]] = r
        logger.info("Resuming: %d experiments already done.", len(completed))

    # Select experiments to run
    experiments = [
        e for e in EXPERIMENTS
        if (args.only is None or e.name in args.only)
        and (not args.resume or e.name not in completed)
    ]
    logger.info("Will run %d experiments.", len(experiments))

    # Load data once
    logger.info("Loading data...")
    final_data, features = load_and_validate_final_data(version="oligo")
    all_val_df = final_data[final_data["split"] == "val"].reset_index(drop=True)
    logger.info("Data loaded | %d train rows | %d val rows | %d features",
                (final_data["split"] == "train").sum(), len(all_val_df), len(features))

    results = list(completed.values())

    for i, cfg in enumerate(experiments):
        logger.info("\n[%d/%d] Starting %s", i + 1, len(experiments), cfg.name)
        try:
            result = run_experiment(cfg, final_data, features, all_val_df,
                                    nthread=args.cpus, seed=args.seed)
            results.append(result)
        except Exception as e:
            logger.error("FAILED: %s — %s", cfg.name, e, exc_info=True)
            results.append({"name": cfg.name, "error": str(e)})

        # Save after every experiment
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)

    # Summary table
    logger.info("\n\n" + "="*70)
    logger.info("SUMMARY (sorted by val Spearman)")
    logger.info("="*70)
    logger.info("%-22s  %8s  %8s  %8s", "experiment", "spearman", "features", "time(s)")
    logger.info("-"*70)
    ranked = sorted(
        [r for r in results if "val_spearman" in r],
        key=lambda r: r["val_spearman"],
        reverse=True,
    )
    for r in ranked:
        logger.info("%-22s  %8.4f  %8d  %8.0f",
                    r["name"], r["val_spearman"], r.get("n_features", -1), r.get("elapsed_s", 0))

    logger.info("\nFull results saved -> %s", args.output)


if __name__ == "__main__":
    main()
