import gc
import json
import logging
import os

import numpy as np
import numpy as cp  # swap for cupy on GPU: import cupy as cp
import optuna
import pandas as pd
import xgboost as xgb
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error
from tqdm.auto import tqdm

from tauso.data.consts import INHIBITION

logger = logging.getLogger(__name__)

PARSIMONY_TOLERANCE = 0.005
EARLY_STOP_DROP = 0.15


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def split_data(df, features, split_col="split", train_val="train", val_val="val", test_val="test"):
    return (
        df[df[split_col] == train_val].copy(),
        df[df[split_col] == val_val].copy(),
        df[df[split_col] == test_val].copy(),
    )


def get_large_cohort_indices(df, group_cols, min_size=50):
    df = df.reset_index(drop=True)
    return [g.index.values for _, g in df.groupby(group_cols) if len(g) >= min_size]


def _collect_group_stats(y_true, preds, eval_groups):
    """Returns (spearmans, top1_arrays, top5_arrays) over eval groups."""
    spearmans, top1s, top5s = [], [], []
    for idxs in eval_groups:
        t, p = y_true[idxs], preds[idxs]
        corr, _ = spearmanr(t, p)
        if not np.isnan(corr):
            spearmans.append(corr)
        n = len(t)
        k1, k5 = max(1, int(n * 0.01)), max(1, int(n * 0.05))
        top1s.append(t[np.argpartition(p, -k1)[-k1:]])
        top5s.append(t[np.argpartition(p, -k5)[-k5:]])
    return spearmans, top1s, top5s


def calculate_metrics(preds, y_true, eval_groups):
    spearmans, top1s, top5s = _collect_group_stats(y_true, preds, eval_groups)
    return {
        "spearman": np.nanmedian(spearmans) if spearmans else 0.0,
        "top1_inhibition": float(np.nanmean([np.mean(t) for t in top1s])) if top1s else 0.0,
        "top5_inhibition": float(np.nanmean([np.mean(t) for t in top5s])) if top5s else 0.0,
    }


def evaluate_final_performance(model, X_np, y_true, eval_groups, feature_names):
    preds = model.predict(xgb.DMatrix(X_np, feature_names=feature_names))
    spearmans, top1s, top5s = _collect_group_stats(y_true, preds, eval_groups)
    return {
        "MAE": float(mean_absolute_error(y_true, preds)),
        "RMSE": float(np.sqrt(mean_squared_error(y_true, preds))),
        "Spearman": float(np.nanmedian(spearmans)) if spearmans else 0.0,
        "Median_Top1_Inhib": float(np.nanmedian([np.median(t) for t in top1s])) if top1s else 0.0,
        "Median_Top5_Inhib": float(np.nanmedian([np.median(t) for t in top5s])) if top5s else 0.0,
    }


def _group_spearman(preds, eval_groups):
    scores = [spearmanr(y, preds[idxs])[0] for idxs, y in eval_groups]
    scores = [s for s in scores if not np.isnan(s)]
    return np.nanmean(scores) if scores else 0.0


def _permutation_importance(bst, X_val, y_val, val_select_idx, current_features, n_repeats=3, seed=0):
    """
    Computes importance as mean Spearman drop on val when each feature is shuffled.
    Higher = more important. Used as a drop-in replacement for gain importance in RFE.
    """
    rng = np.random.default_rng(seed)
    baseline = calculate_metrics(bst.inplace_predict(X_val), y_val, val_select_idx)["spearman"]
    importances = {}
    for i, feat in enumerate(current_features):
        drops = []
        for _ in range(n_repeats):
            X_perturbed = X_val.copy()
            X_perturbed[:, i] = rng.permutation(X_perturbed[:, i])
            preds = bst.inplace_predict(X_perturbed)
            if hasattr(preds, "get"):
                preds = preds.get()
            perturbed_spear = calculate_metrics(preds, y_val, val_select_idx)["spearman"]
            drops.append(baseline - perturbed_spear)
        importances[feat] = float(np.mean(drops))
    return importances


def _log_trial(study, trial):
    logger.info(
        "Optuna trial %d | spearman=%.4f | best=%.4f | rounds=%s depth=%s lr=%.4f",
        trial.number + 1,
        trial.value if trial.value is not None else float("nan"),
        study.best_value,
        trial.params.get("num_boost_round"),
        trial.params.get("max_depth"),
        trial.params.get("learning_rate", 0.0),
    )


def _save_outputs(path_base, model_tv, model_all, metrics_tv, metrics_all, history_df):
    path_tv = path_base.replace(".json", "_TrainVal.json")
    path_ad = path_base.replace(".json", "_AllData.json")
    model_tv.save_model(path_tv)
    model_all.save_model(path_ad)
    for path, metrics in [(path_tv, metrics_tv), (path_ad, metrics_all)]:
        with open(path.replace(".json", "_metrics.json"), "w") as f:
            json.dump(metrics, f, indent=4)
    history_df.to_csv(path_base.replace(".json", "_RFE_History.csv"), index=False)
    logger.info("TrainVal -> %s", path_tv)
    logger.info("AllData  -> %s", path_ad)


# ---------------------------------------------------------------------------
# Training steps
# ---------------------------------------------------------------------------

def tune_hyperparameters(train_data, val_data, objective, val_eval_groups_optuna, seed, device="cpu", n_jobs=1):
    X_train, y_train = train_data
    X_val, y_val = val_data
    dtrain = xgb.DMatrix(X_train, label=y_train)
    dval = xgb.DMatrix(X_val, label=y_val)

    def objective_func(trial):
        params = {
            "tree_method": "hist", "device": device, "nthread": n_jobs,
            "objective": objective, "seed": seed,
            "max_depth": trial.suggest_int("max_depth", 2, 6),
            "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.1, log=True),
            "subsample": trial.suggest_float("subsample", 0.4, 0.9),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.3, 0.9),
            "min_child_weight": trial.suggest_int("min_child_weight", 50, 300),
            "gamma": trial.suggest_float("gamma", 0.1, 10.0, log=True),
            "reg_alpha": trial.suggest_float("reg_alpha", 0.1, 100.0, log=True),
            "reg_lambda": trial.suggest_float("reg_lambda", 0.1, 100.0, log=True),
        }
        if objective == "reg:pseudohubererror":
            params["huber_slope"] = trial.suggest_float("huber_slope", 0.5, 10.0, log=True)
        bst = xgb.train(
            params, dtrain,
            num_boost_round=trial.suggest_int("num_boost_round", 200, 1500, step=100),
            evals=[(dval, "val")], early_stopping_rounds=50, verbose_eval=False,
        )
        trial.set_user_attr("best_iteration", bst.best_iteration)
        return _group_spearman(bst.predict(dval), val_eval_groups_optuna)

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=seed))
    logger.info("Optuna tuning: 50 trials | objective=%s device=%s n_jobs=%d", objective, device, n_jobs)
    study.optimize(objective_func, n_trials=50, show_progress_bar=True, callbacks=[_log_trial])

    best_params = study.best_params.copy()
    best_params["num_boost_round"] = study.best_trial.user_attrs["best_iteration"]
    best_params.update({"tree_method": "hist", "device": device, "nthread": n_jobs, "objective": objective, "seed": seed})

    logger.info("Tuning complete | best_spearman=%.4f | params=%s", study.best_value, best_params)
    return best_params


def run_backward_selection(
    train_data, val_data, test_data, features, best_params, val_eval_idx, val_select_idx, test_select_idx, seed,
    importance_type="gain", parsimony_tolerance=PARSIMONY_TOLERANCE,
):
    X_train, y_train, y_train_np = train_data
    X_val, y_val = val_data
    X_test, y_test = test_data

    feat_to_idx = {f: i for i, f in enumerate(features)}
    current_features = list(features)
    history, max_spearman = [], -1.0

    params = best_params.copy()
    if params.get("seed") != seed:
        raise ValueError(f"Seed mismatch: best_params={params.get('seed')} vs argument={seed}")
    num_rounds = params.pop("num_boost_round", 1000)

    logger.info("RFE start | features=%d rounds=%d importance=%s parsimony=%.4f",
                len(current_features), num_rounds, importance_type, parsimony_tolerance)
    total = ((len(current_features) - 100) // 5) + 100 if len(current_features) > 100 else len(current_features)
    pbar = tqdm(total=total, desc="Dropping Features")

    for step in range(len(current_features)):
        idxs = [feat_to_idx[f] for f in current_features]
        Xtr, Xvl, Xte = X_train[:, idxs], X_val[:, idxs], X_test[:, idxs]

        dtrain = xgb.QuantileDMatrix(Xtr, label=y_train, feature_names=current_features)
        bst = xgb.train(params, dtrain, num_boost_round=num_rounds, verbose_eval=False)

        tr_preds, vl_preds, te_preds = bst.inplace_predict(Xtr), bst.inplace_predict(Xvl), bst.inplace_predict(Xte)
        if hasattr(tr_preds, "get"):
            tr_preds, vl_preds, te_preds = tr_preds.get(), vl_preds.get(), te_preds.get()

        m_train = calculate_metrics(tr_preds, y_train_np, val_eval_idx)
        val_spear = calculate_metrics(vl_preds, y_val, val_select_idx)["spearman"]
        max_spearman = max(max_spearman, val_spear)

        logger.info(
            "RFE step %d | features=%d | val=%.4f train=%.4f peak=%.4f delta=%.4f",
            step, len(current_features), val_spear, m_train["spearman"], max_spearman, val_spear - max_spearman,
        )

        row = {
            "num_features": len(current_features),
            "features_list": current_features.copy(),
            "train_spearman": m_train["spearman"],
            "val_spearman": calculate_metrics(vl_preds, y_val, val_eval_idx)["spearman"],
            "val_spearman_select": val_spear,
        }

        if val_spear < max_spearman - EARLY_STOP_DROP:
            logger.info("Early stop: val %.4f dropped %.2f below peak %.4f", val_spear, EARLY_STOP_DROP, max_spearman)
            row["dropped_features"] = []
            history.append(row)
            break

        if len(current_features) == 1:
            row["dropped_features"] = current_features
            history.append(row)
            pbar.update(1)
            break

        if importance_type == "permutation":
            imp = _permutation_importance(bst, Xvl, y_val, val_select_idx, current_features, seed=seed)
        else:
            raw = bst.get_score(importance_type="gain")
            imp = {f: raw.get(f, 0.0) for f in current_features}
        ranked = sorted(imp.items(), key=lambda x: x[1])
        drop_n = min(5 if len(current_features) > 100 else 1, len(current_features) - 1)
        to_drop = [f for f, _ in ranked[:drop_n]]
        logger.debug("Dropping %d: %s", drop_n, to_drop)

        row["dropped_features"] = to_drop
        history.append(row)
        current_features = [f for f in current_features if f not in to_drop]
        pbar.update(1)

        del dtrain, bst, Xtr, Xvl, Xte
        if hasattr(cp, "get_default_memory_pool"):
            cp.get_default_memory_pool().free_all_blocks()
        gc.collect()

    pbar.close()

    df_hist = pd.DataFrame(history)
    threshold = df_hist["val_spearman_select"].max() - parsimony_tolerance
    optimal = df_hist[df_hist["val_spearman_select"] >= threshold].iloc[-1]["features_list"]

    logger.info("RFE complete | selected=%d / %d | best_val=%.4f", len(optimal), len(features), df_hist["val_spearman_select"].max())
    return optimal, df_hist


def train_final_model(optimal_features, features, train_data, best_params, seed):
    X, y = train_data
    feat_to_idx = {f: i for i, f in enumerate(features)}
    idxs = [feat_to_idx[f] for f in optimal_features]
    params = {**best_params, "seed": seed}
    num_rounds = params.pop("num_boost_round", 1000)

    logger.info("Training final model | features=%d rounds=%d", len(optimal_features), num_rounds)
    dtrain = xgb.QuantileDMatrix(X[:, idxs], label=y, feature_names=optimal_features)
    return xgb.train(params, dtrain, num_boost_round=num_rounds, verbose_eval=False)


# ---------------------------------------------------------------------------
# Convenience orchestrator (used by notebooks)
# ---------------------------------------------------------------------------

def run_pipeline(
    final_data, features, split_config, eval_group, select_group,
    loss_type="L1", save_path="model.json", seed=1, device="cpu", n_jobs=1,
    importance_type="gain", parsimony_tolerance=PARSIMONY_TOLERANCE,
):
    os.environ["PYTHONHASHSEED"] = str(seed)
    _OBJECTIVES = {"L1": "reg:absoluteerror", "L2": "reg:squarederror", "huber": "reg:pseudohubererror"}
    if loss_type not in _OBJECTIVES:
        raise ValueError(f"loss_type must be one of {list(_OBJECTIVES)}. Got: {loss_type!r}")
    objective = _OBJECTIVES[loss_type]
    logger.info("Pipeline | objective=%s eval=%s seed=%d", objective, eval_group, seed)

    final_data = final_data.sort_values(by=features).reset_index(drop=True)
    train_df, val_df, test_df = split_data(final_data, features, split_col=split_config["col"],
                                           train_val=split_config["train"], val_val=split_config["val"],
                                           test_val=split_config["test"])

    X_train, y_train = cp.array(train_df[features].values), cp.array(train_df[INHIBITION].values)
    X_val, y_val     = cp.array(val_df[features].values),   val_df[INHIBITION].values
    X_test, y_test   = cp.array(test_df[features].values),  test_df[INHIBITION].values
    y_train_np       = train_df[INHIBITION].values

    train_eval_idx  = get_large_cohort_indices(train_df, eval_group)
    val_eval_idx    = get_large_cohort_indices(val_df, eval_group)
    test_eval_idx   = get_large_cohort_indices(test_df, eval_group)
    val_select_idx  = get_large_cohort_indices(val_df, select_group)
    test_select_idx = get_large_cohort_indices(test_df, select_group)

    val_eval_groups_optuna = [
        (g.index.values, y_val[g.index.values])
        for _, g in val_df.reset_index(drop=True).groupby(eval_group)
        if len(g) >= 20
    ]

    logger.info("[1/3] Hyperparameter tuning...")
    best_params = tune_hyperparameters(
        (train_df[features].values, y_train_np), (val_df[features].values, y_val),
        objective, val_eval_groups_optuna, seed, device=device, n_jobs=n_jobs,
    )

    logger.info("[2/3] Backward feature selection...")
    optimal_features, history_df = run_backward_selection(
        (X_train, y_train, y_train_np), (X_val, y_val), (X_test, y_test),
        features, best_params, val_eval_idx, val_select_idx, test_select_idx, seed,
        importance_type=importance_type, parsimony_tolerance=parsimony_tolerance,
    )

    logger.info("[3/3] Training final models...")
    X_tv = cp.vstack((X_train, X_val))
    y_tv = cp.concatenate((y_train, cp.array(y_val)))
    X_all = cp.vstack((X_tv, X_test))
    y_all = cp.concatenate((y_tv, cp.array(y_test)))

    model_tv  = train_final_model(optimal_features, features, (X_tv, y_tv), best_params, seed)
    model_all = train_final_model(optimal_features, features, (X_all, y_all), best_params, seed)

    feat_to_idx = {f: i for i, f in enumerate(features)}
    opt_idxs = [feat_to_idx[f] for f in optimal_features]
    X_tr_opt, X_vl_opt = train_df[features].values[:, opt_idxs], val_df[features].values[:, opt_idxs]

    metrics_tv = {
        "Train":      evaluate_final_performance(model_tv, X_tr_opt, y_train_np, train_eval_idx, optimal_features),
        "Validation": evaluate_final_performance(model_tv, X_vl_opt, y_val,      val_eval_idx,   optimal_features),
    }
    metrics_all = {
        "Train":      evaluate_final_performance(model_all, X_tr_opt, y_train_np, train_eval_idx, optimal_features),
        "Validation": evaluate_final_performance(model_all, X_vl_opt, y_val,      val_eval_idx,   optimal_features),
    }

    _save_outputs(save_path, model_tv, model_all, metrics_tv, metrics_all, history_df)
    logger.info("Pipeline complete.")
    return {"model_train_val": model_tv, "model_all": model_all}, optimal_features, history_df
