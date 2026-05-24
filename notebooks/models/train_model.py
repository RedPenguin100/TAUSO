"""
Train the XGBoost ASO efficacy model on server CPUs (or GPU if available).

Steps can be run together or independently so that long runs can be split
across multiple jobs or resumed after a failure:

  tune   -- Optuna hyperparameter search. Saves best_params_<name>.json.
  select -- Backward feature selection (RFE). Loads best_params, saves
             optimal_features_<name>.json and RFE history CSV.
  train  -- Train final models on selected features. Loads best_params and
             optimal_features, saves TrainVal / AllData model JSONs + metrics.
  all    -- Runs tune -> select -> train in sequence (default).

Usage:
    python -m notebooks.models.train_model --loss L2 --split cohort --cpus 32
    python -m notebooks.models.train_model --loss L2 --split cohort --step tune --cpus 32
    python -m notebooks.models.train_model --loss L2 --split cohort --step select --cpus 32
    python -m notebooks.models.train_model --loss L2 --split cohort --step train
    python -m notebooks.models.train_model --loss L2 --split cohort --device cuda
    python -m notebooks.models.train_model --loss L2 --split cohort --cpus 32 --seed 42

Output is written to notebooks/models/SeenOligoModel/ by default (--output to override).

SLURM note: run with `python -u` or PYTHONUNBUFFERED=1 for live log output.
"""
import argparse
import json
import logging
import os
import re
import sys
from pathlib import Path

import numpy as np
from notebooks.models.utility import load_and_validate_final_data

from tauso.data.consts import CANONICAL_GENE, CELL_LINE, INHIBITION

logger = logging.getLogger(__name__)

DEFAULT_OUTPUT_DIR = Path(__file__).parent / "SeenOligoModel"

SPLIT_CONFIGS = {
    "cohort":    {"eval_group": "custom_id", "select_group": [CANONICAL_GENE, CELL_LINE], "suffix": ""},
    "custom_id": {"eval_group": "custom_id", "select_group": ["custom_id"],               "suffix": "_CustomId"},
}

SPLIT_CONFIG = {"col": "split", "train": "train", "val": "val", "test": "test"}


# ---------------------------------------------------------------------------
# File I/O helpers
# ---------------------------------------------------------------------------

def _load_json(path, missing_step):
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run --step {missing_step} first.")
    with open(path) as f:
        return json.load(f)


def _save_json(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=4)


def _apply_feature_filter(features, pattern):
    """Exclude features whose name matches the regex pattern."""
    if not pattern:
        return features
    rx = re.compile(pattern)
    keep    = [f for f in features if not rx.search(f)]
    dropped = [f for f in features if     rx.search(f)]
    if dropped:
        logger.info(f"Feature filter | excluding {len(dropped)} features matching '{pattern}' (examples: {dropped[:5]})")
    return keep


def _drop_zero_variance(final_data, features, dropped_log_path=None):
    """Drop features with <=1 unique non-NaN value (cannot inform tree splits).

    Returns the kept feature list. If dropped_log_path is given and anything was
    dropped, writes a JSON mapping {feature: constant_value} ("NaN" for all-NaN
    columns) so the dropped set is inspectable after the run.
    """
    dropped = {}
    keep = []
    for f in features:
        unique = final_data[f].dropna().unique()
        if len(unique) > 1:
            keep.append(f)
        elif len(unique) == 1:
            dropped[f] = float(unique[0])
        else:
            dropped[f] = "NaN"

    if dropped:
        logger.info("Zero-variance filter: dropped %d / %d features: %s",
                    len(dropped), len(features), list(dropped.keys()))
        if dropped_log_path is not None:
            _save_json(dropped, dropped_log_path)
            logger.info("Dropped-features log -> %s", dropped_log_path)
    return keep


def _resolve_features(final_data, all_features, paths):
    """Return filtered feature list, saving to disk on first call and loading on subsequent ones.

    If the cache exists but is stale (columns missing from data, or new columns in data
    that weren't there when the cache was written), the cache is automatically regenerated.
    The zero-variance filter is applied on every call so newly-constant columns get
    dropped even when the cache is reused.
    """
    if paths["input_features"].exists():
        cached = _load_json(paths["input_features"], "tune")
        missing_from_data = [f for f in cached if f not in all_features]
        new_in_data = [f for f in all_features if f not in cached]

        if missing_from_data or new_in_data:
            logger.warning(
                "Cached input_features is stale: %d columns gone from data, "
                "%d new columns in data. Regenerating %s.",
                len(missing_from_data), len(new_in_data), paths["input_features"],
            )
            paths["input_features"].unlink()
        else:
            logger.info("Loaded %d input features from %s", len(cached), paths["input_features"])
            features = _drop_zero_variance(final_data, cached, paths["dropped_features"])
            if len(features) != len(cached):
                _save_json(features, paths["input_features"])
                logger.info("Cache updated after zero-variance filter -> %s", paths["input_features"])
            return features

    features = _drop_zero_variance(final_data, all_features, paths["dropped_features"])
    _save_json(features, paths["input_features"])
    logger.info("%d input features saved -> %s", len(features), paths["input_features"])
    return features


def _save_model_and_metrics(model, path, metrics):
    model.save_model(path)
    _save_json(metrics, path.replace(".json", "_metrics.json"))
    logger.info("  %s | train_spearman=%.4f val_spearman=%.4f",
                Path(path).name, metrics["Train"]["Spearman"], metrics["Validation"]["Spearman"])


def _intermediate_paths(args):
    name = f"{args.loss}{SPLIT_CONFIGS[args.split]['suffix']}"
    return {
        "name":             name,
        "model_base":       str(args.output / f"Model_Oligo_{name}.json"),
        "input_features":   args.output / f"input_features_{name}.json",
        "dropped_features": args.output / f"dropped_features_{name}.json",
        "best_params":      args.output / f"best_params_{name}.json",
        "optimal_features": args.output / f"optimal_features_{name}.json",
        "rfe_history":      args.output / f"Model_Oligo_{name}_RFE_History.csv",
    }


# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def _prep_data(final_data, features, eval_group, select_group, min_eval_size=50, rank=False):
    from notebooks.models.SeenOligoModel.base_model import get_large_cohort_indices, split_data

    # Learning-to-rank needs rows grouped (contiguous) by query id. Sort the whole
    # frame by the ranking cohort (select_group) so each split stays grouped, and
    # assign a global qid (consistent across splits so stacking can re-group later).
    if rank:
        final_data = final_data.sort_values(by=list(select_group)).reset_index(drop=True)
        qid_all = final_data.groupby(select_group, sort=False).ngroup().values.astype(np.int32)
        final_data = final_data.assign(_qid=qid_all)
    else:
        final_data = final_data.sort_values(by=features).reset_index(drop=True)

    train_df, val_df, test_df = split_data(
        final_data, features,
        split_col=SPLIT_CONFIG["col"], train_val=SPLIT_CONFIG["train"],
        val_val=SPLIT_CONFIG["val"], test_val=SPLIT_CONFIG["test"],
    )

    def arrays(df):
        return df[features].values.astype(np.float32), df[INHIBITION].values.astype(np.float32)

    X_train, y_train = arrays(train_df)
    X_val,   y_val   = arrays(val_df)
    X_test,  y_test  = arrays(test_df)

    logger.info("Split | train=%d val=%d test=%d | min_eval_size=%d", len(train_df), len(val_df), len(test_df), min_eval_size)

    val_optuna_groups = [
        (g.index.values, y_val[g.index.values])
        for _, g in val_df.reset_index(drop=True).groupby(eval_group)
        if len(g) >= min_eval_size
    ]

    out = {
        "X_train": X_train, "y_train": y_train,
        "X_val":   X_val,   "y_val":   y_val,
        "X_test":  X_test,  "y_test":  y_test,
        "train_eval_idx":       get_large_cohort_indices(train_df, eval_group, min_size=min_eval_size),
        "val_eval_idx":         get_large_cohort_indices(val_df,   eval_group, min_size=min_eval_size),
        "test_eval_idx":        get_large_cohort_indices(test_df,  eval_group, min_size=min_eval_size),
        "val_select_idx":       get_large_cohort_indices(val_df,   select_group, min_size=min_eval_size),
        "test_select_idx":      get_large_cohort_indices(test_df,  select_group, min_size=min_eval_size),
        "val_eval_groups_optuna": val_optuna_groups,
    }
    if rank:
        out["qid_train"] = train_df["_qid"].values.astype(np.int32)
        out["qid_val"]   = val_df["_qid"].values.astype(np.int32)
        out["qid_test"]  = test_df["_qid"].values.astype(np.int32)
    return out


# ---------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------

def step_tune(data, args, paths, objective):
    from notebooks.models.SeenOligoModel.base_model import tune_hyperparameters
    logger.info("=== STEP: tune ===")
    best_params = tune_hyperparameters(
        (data["X_train"], data["y_train"]), (data["X_val"], data["y_val"]),
        objective, data["val_eval_groups_optuna"],
        seed=args.seed, device=args.device, n_jobs=args.cpus,
        qid_train=data.get("qid_train"), qid_val=data.get("qid_val"),
    )
    _save_json(best_params, paths["best_params"])
    logger.info("Best params saved -> %s", paths["best_params"])
    return best_params


def step_select(data, features, args, paths):
    from notebooks.models.SeenOligoModel.base_model import run_backward_selection
    logger.info("=== STEP: select ===")
    features = _load_json(paths["input_features"], "tune")
    best_params = _load_json(paths["best_params"], "tune")

    optimal, history_df = run_backward_selection(
        (data["X_train"], data["y_train"], data["y_train"]),
        (data["X_val"],   data["y_val"]),
        (data["X_test"],  data["y_test"]),
        features, best_params,
        data["val_eval_idx"], data["val_select_idx"], data["test_select_idx"],
        seed=args.seed,
        importance_type=args.importance,
        parsimony_tolerance=args.parsimony_tolerance,
        device=args.device,
        history_path=paths["rfe_history"],
        optimal_path=paths["optimal_features"],
        qid_train=data.get("qid_train"),
    )
    logger.info("Optimal features (%d) -> %s", len(optimal), paths["optimal_features"])
    logger.info("RFE history -> %s", paths["rfe_history"])
    return optimal


def step_train(data, features, args, paths):
    from notebooks.models.SeenOligoModel.base_model import evaluate_final_performance, train_final_model
    logger.info("=== STEP: train ===")

    features         = _load_json(paths["input_features"],   "tune")
    best_params      = _load_json(paths["best_params"],      "tune")
    optimal_features = _load_json(paths["optimal_features"], "select")
    logger.info("Loaded %d features | rounds=%s", len(optimal_features), best_params.get("num_boost_round"))

    is_rank = str(best_params.get("objective", "")).startswith("rank:")

    def stack(*keys):
        X = np.vstack([data[f"X_{k}"] for k in keys])
        y = np.concatenate([data[f"y_{k}"] for k in keys])
        if not is_rank:
            return X, y, None
        # ranking needs contiguous qid groups; concatenated splits interleave the
        # same cohort, so re-sort the stacked rows by qid (stable keeps determinism)
        q = np.concatenate([data[f"qid_{k}"] for k in keys])
        order = np.argsort(q, kind="stable")
        return X[order], y[order], q[order]

    X_tv,  y_tv,  q_tv  = stack("train", "val")
    X_all, y_all, q_all = stack("train", "val", "test")

    def train(X, y, q, label):
        logger.info("Training %s (%d rows)", label, len(X))
        return train_final_model(optimal_features, features, (X, y), best_params, seed=args.seed, qid=q)

    model_tv  = train(X_tv,  y_tv,  q_tv,  "train+val")
    model_all = train(X_all, y_all, q_all, "all data")

    opt_idxs = [features.index(f) for f in optimal_features]

    def eval_split(model, X_key, y_key, idx_key):
        return evaluate_final_performance(
            model, data[X_key][:, opt_idxs], data[y_key], data[idx_key], optimal_features
        )

    for model, suffix, _label in [(model_tv, "_TrainVal", "TrainVal"), (model_all, "_AllData", "AllData")]:
        metrics = {
            "Train":      eval_split(model, "X_train", "y_train", "train_eval_idx"),
            "Validation": eval_split(model, "X_val",   "y_val",   "val_eval_idx"),
        }
        path = paths["model_base"].replace(".json", f"{suffix}.json")
        _save_model_and_metrics(model, path, metrics)

    return {"model_train_val": model_tv, "model_all": model_all}


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--loss",      required=True, choices=["L1", "L2", "huber", "ndcg"])
    parser.add_argument("--split",     required=True, choices=list(SPLIT_CONFIGS))
    parser.add_argument("--step",      default="all", choices=["all", "tune", "select", "train"])
    parser.add_argument("--cpus",      type=int,  default=1)
    parser.add_argument("--device",    default="cpu", choices=["cpu", "cuda"])
    parser.add_argument("--seed",      type=int,  default=1)
    parser.add_argument("--output",    type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--importance", default="gain", choices=["gain", "permutation"],
                        help="Feature importance metric for RFE (default: gain). "
                             "permutation is slower but more reliable.")
    parser.add_argument("--parsimony-tolerance", type=float, default=0.005,
                        help="Max Spearman drop from peak allowed when picking parsimonious feature set "
                             "(default: 0.005).")
    parser.add_argument("--split-source", default="oligoai", choices=["oligoai", "tauso"],
                        help="'oligoai' uses the existing split column (default); "
                             "'tauso' creates a stratified temporal split per gene×cell-line.")
    parser.add_argument("--min-eval-size", type=int, default=50,
                        help="Minimum cohort size to include in evaluation metrics (default: 50). "
                             "Floored at 20 — smaller cohorts give noisy Spearman estimates.")
    parser.add_argument("--exclude-feature-pattern", default=None,
                        help="Regex; features whose name matches are excluded from training. "
                             "Example: --exclude-feature-pattern '^OHE_pos' drops all positional one-hots.")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s | %(message)s",
        stream=sys.stdout, force=True,
    )

    if args.min_eval_size < 20:
        logger.warning("--min-eval-size %d is below floor of 20; clamping to 20.", args.min_eval_size)
        args.min_eval_size = 20

    os.environ["PYTHONHASHSEED"] = str(args.seed)
    args.output.mkdir(parents=True, exist_ok=True)

    paths    = _intermediate_paths(args)
    cfg      = SPLIT_CONFIGS[args.split]
    _OBJECTIVES = {"L1": "reg:absoluteerror", "L2": "reg:squarederror",
                   "huber": "reg:pseudohubererror", "ndcg": "rank:ndcg"}
    objective = _OBJECTIVES[args.loss]
    is_rank = objective.startswith("rank:")

    logger.info("loss=%s split=%s step=%s device=%s cpus=%d seed=%d", args.loss, args.split, args.step, args.device, args.cpus, args.seed)

    logger.info("Loading data...")
    final_data, all_features = load_and_validate_final_data(version="oligo", split_source=args.split_source)
    features = _resolve_features(final_data, all_features, paths)
    if args.exclude_feature_pattern:
        before = len(features)
        features = _apply_feature_filter(features, args.exclude_feature_pattern)
        if len(features) != before:
            _save_json(features, paths["input_features"])
    logger.info("Features: %d", len(features))

    data = _prep_data(final_data, features, cfg["eval_group"], cfg["select_group"],
                      min_eval_size=args.min_eval_size, rank=is_rank)

    if args.step in ("tune",   "all"): step_tune(data, args, paths, objective)
    if args.step in ("select", "all"): step_select(data, features, args, paths)
    if args.step in ("train",  "all"): step_train(data, features, args, paths)

    logger.info("Done.")


if __name__ == "__main__":
    main()
