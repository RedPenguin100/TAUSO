import gc
import json
import os

import pandas as pd
import numpy as np
import cupy as cp
import xgboost as xgb
import optuna
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error
from tqdm.auto import tqdm

from tauso.data.consts import INHIBITION

# Constants
PARSIMONY_TOLERANCE = 0.03
EARLY_STOP_DROP = 0.15


def split_data(df, features, split_col="split", train_val="train", val_val="val", test_val="test"):
    train_df = df[df[split_col] == train_val].copy()
    val_df = df[df[split_col] == val_val].copy()
    test_df = df[df[split_col] == test_val].copy()
    return train_df, val_df, test_df


def get_large_cohort_indices(df, group_cols, min_size=50):
    indices = []
    df_reset = df.reset_index(drop=True)
    for cohort, group in df_reset.groupby(group_cols):
        if len(group) >= min_size:
            indices.append(group.index.values)
    return indices


def calculate_metrics(preds, y_true, eval_groups):
    spearmans, top_1_means, top_5_means = [], [], []
    for idxs in eval_groups:
        t_vals, p_vals = y_true[idxs], preds[idxs]
        corr, _ = spearmanr(t_vals, p_vals)
        if not np.isnan(corr):
            spearmans.append(corr)

        n = len(t_vals)
        k1, k5 = max(1, int(n * 0.01)), max(1, int(n * 0.05))

        if k5 > 0:
            top_5_means.append(np.mean(t_vals[np.argpartition(p_vals, -k5)[-k5:]]))
        if k1 > 0:
            top_1_means.append(np.mean(t_vals[np.argpartition(p_vals, -k1)[-k1:]]))

    return {
        "spearman": np.nanmedian(spearmans) if spearmans else 0.0,
        "top1_inhibition": np.nanmean(top_1_means) if top_1_means else 0.0,
        "top5_inhibition": np.nanmean(top_5_means) if top_5_means else 0.0,
    }


def evaluate_final_performance(model, X_np, y_true, eval_groups, feature_names):
    """Calculates specific final metrics for JSON tracking without altering RFE tracking."""
    dmat = xgb.DMatrix(X_np, feature_names=feature_names)
    preds = model.predict(dmat)

    mae = mean_absolute_error(y_true, preds)
    rmse = np.sqrt(mean_squared_error(y_true, preds))

    spearmans, top_1_medians, top_5_medians = [], [], []
    for idxs in eval_groups:
        t_vals, p_vals = y_true[idxs], preds[idxs]

        corr, _ = spearmanr(t_vals, p_vals)
        if not np.isnan(corr):
            spearmans.append(corr)

        n = len(t_vals)
        k1, k5 = max(1, int(n * 0.01)), max(1, int(n * 0.05))

        if k5 > 0:
            top_5_medians.append(np.median(t_vals[np.argpartition(p_vals, -k5)[-k5:]]))
        if k1 > 0:
            top_1_medians.append(np.median(t_vals[np.argpartition(p_vals, -k1)[-k1:]]))

    return {
        "MAE": float(mae),
        "RMSE": float(rmse),
        "Spearman": float(np.nanmedian(spearmans)) if spearmans else 0.0,
        "Median_Top1_Inhib": float(np.nanmedian(top_1_medians)) if top_1_medians else 0.0,
        "Median_Top5_Inhib": float(np.nanmedian(top_5_medians)) if top_5_medians else 0.0,
    }


def tune_hyperparameters(train_data, val_data, objective, val_eval_groups_optuna, seed):
    X_train_np, y_train_np = train_data
    X_val_np, y_val_np = val_data

    dtrain_opt = xgb.DMatrix(X_train_np, label=y_train_np)
    dval_opt = xgb.DMatrix(X_val_np, label=y_val_np)

    def calculate_fast_spearman(preds, eval_groups):
        spearmans = [
            spearmanr(true_vals, preds[idxs])[0]
            for idxs, true_vals in eval_groups
            if not np.isnan(spearmanr(true_vals, preds[idxs])[0])
        ]
        return np.nanmean(spearmans) if spearmans else 0.0

    def objective_func(trial):
        params = {
            "tree_method": "hist",
            "device": "cuda",
            "objective": objective,
            "seed": seed,  # <-- Local seed for XGBoost
            "max_depth": trial.suggest_int("max_depth", 2, 6),
            "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.1, log=True),
            "subsample": trial.suggest_float("subsample", 0.4, 0.9),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.3, 0.9),
            "min_child_weight": trial.suggest_int("min_child_weight", 50, 300),
            "gamma": trial.suggest_float("gamma", 0.1, 10.0, log=True),
            "reg_alpha": trial.suggest_float("reg_alpha", 0.1, 100.0, log=True),
            "reg_lambda": trial.suggest_float("reg_lambda", 0.1, 100.0, log=True),
        }

        bst = xgb.train(
            params,
            dtrain_opt,
            num_boost_round=trial.suggest_int("num_boost_round", 200, 1500, step=100),
            evals=[(dval_opt, "val")],
            early_stopping_rounds=50,
            verbose_eval=False,
        )

        # <-- FIX 2: Save the ACTUAL best iteration where early stopping halted
        trial.set_user_attr("best_iteration", bst.best_iteration)

        return calculate_fast_spearman(bst.predict(dval_opt), val_eval_groups_optuna)

    sampler = optuna.samplers.TPESampler(seed=seed)  # <-- Local seed for Optuna
    study = optuna.create_study(direction="maximize", sampler=sampler)
    study.optimize(objective_func, n_trials=50, show_progress_bar=True)

    best_params = study.best_params.copy()

    actual_best_rounds = study.best_trial.user_attrs["best_iteration"]
    best_params["num_boost_round"] = actual_best_rounds

    best_params.update(
        {"tree_method": "hist", "device": "cuda", "objective": objective, "seed": seed}
    )

    return best_params


# === SEED ADDED AS PARAMETER ===
def run_backward_selection(
    train_data, val_data, test_data, features, best_params, val_eval_idx, val_select_idx, test_select_idx, seed
):
    X_train_gpu, y_train_gpu, y_train_np = train_data
    X_val_gpu, y_val_np = val_data
    X_test_gpu, y_test_np = test_data

    feat_to_idx = {f: i for i, f in enumerate(features)}
    current_features = list(features)
    selection_history = []
    max_spearman_seen = -1.0

    selection_params = best_params.copy()

    # 1. Ensure "seed" exists in the dictionary and matches the expected parameter
    if "seed" not in selection_params:
        raise KeyError("Critical Error: 'seed' was not found in best_params. Reproducibility is compromised.")

    # 2. Optional: Verify the value matches the function argument to prevent accidental mismatches
    if selection_params["seed"] != seed:
        raise ValueError(f"Mismatched Seed: best_params has {selection_params['seed']} but function received {seed}.")

    # 3. Extract the boost rounds
    num_rounds = selection_params.pop("num_boost_round", 1000)

    total_iters = ((len(current_features) - 100) // 5) + 100 if len(current_features) > 100 else len(current_features)
    pbar = tqdm(total=total_iters, desc="Dropping Features")

    while len(current_features) > 0:
        curr_idxs = [feat_to_idx[f] for f in current_features]
        curr_X_train = X_train_gpu[:, curr_idxs]
        curr_X_val = X_val_gpu[:, curr_idxs]
        curr_X_test = X_test_gpu[:, curr_idxs]

        dtrain = xgb.QuantileDMatrix(curr_X_train, label=y_train_gpu, feature_names=current_features)
        bst = xgb.train(selection_params, dtrain, num_boost_round=num_rounds, verbose_eval=False)

        train_preds, val_preds, test_preds = (
            bst.inplace_predict(curr_X_train),
            bst.inplace_predict(curr_X_val),
            bst.inplace_predict(curr_X_test),
        )
        if hasattr(train_preds, "get"):
            train_preds, val_preds, test_preds = train_preds.get(), val_preds.get(), test_preds.get()

        metrics_train = calculate_metrics(train_preds, y_train_np, val_eval_idx)
        metrics_val = calculate_metrics(val_preds, y_val_np, val_eval_idx)

        val_spear_select = calculate_metrics(val_preds, y_val_np, val_select_idx)["spearman"]

        if val_spear_select > max_spearman_seen:
            max_spearman_seen = val_spear_select

        step_data = {
            "num_features": len(current_features),
            "features_list": current_features.copy(),
            "train_spearman": metrics_train["spearman"],
            "val_spearman": metrics_val["spearman"],
            "val_spearman_select": val_spear_select,
        }

        if val_spear_select < (max_spearman_seen - EARLY_STOP_DROP):
            step_data["dropped_features"] = []
            selection_history.append(step_data)
            break

        if len(current_features) == 1:
            step_data["dropped_features"] = current_features
            selection_history.append(step_data)
            pbar.update(1)
            break

        imp_dict = bst.get_score(importance_type="gain")
        sorted_feats = sorted({f: imp_dict.get(f, 0.0) for f in current_features}.items(), key=lambda x: x[1])

        drop_n = min(5 if len(current_features) > 100 else 1, len(current_features) - 1)
        feats_to_drop = [x[0] for x in sorted_feats[:drop_n]]

        step_data["dropped_features"] = feats_to_drop
        selection_history.append(step_data)
        current_features = [f for f in current_features if f not in feats_to_drop]
        pbar.update(1)

        del dtrain, bst, curr_X_train, curr_X_val, curr_X_test
        cp.get_default_memory_pool().free_all_blocks()
        gc.collect()

    pbar.close()

    df_backward_selection = pd.DataFrame(selection_history)
    threshold = df_backward_selection["val_spearman_select"].max() - PARSIMONY_TOLERANCE
    best_row = df_backward_selection[df_backward_selection["val_spearman_select"] >= threshold].iloc[-1]

    return best_row["features_list"], df_backward_selection


# === SEED ADDED AS PARAMETER ===
def train_final_model(optimal_features, features, train_data, best_params, seed):
    X_gpu, y_gpu = train_data

    feat_to_idx = {f: i for i, f in enumerate(features)}
    opt_idxs = [feat_to_idx[f] for f in optimal_features]

    final_X = X_gpu[:, opt_idxs]
    final_dtrain = xgb.QuantileDMatrix(final_X, label=y_gpu, feature_names=optimal_features)

    train_params = best_params.copy()
    num_rounds = train_params.pop("num_boost_round", 1000)
    train_params["seed"] = seed  # <-- Enforce seed explicitly here

    final_bst = xgb.train(train_params, final_dtrain, num_boost_round=num_rounds, verbose_eval=False)
    return final_bst


def run_pipeline(
    final_data, features, split_config, eval_group, select_group, loss_type="L1", save_path="model.json", seed=1
):
    # <-- FIX 4: Python's hash seed must be set before any sets/dicts are iterated over
    os.environ["PYTHONHASHSEED"] = str(seed)

    # # <-- FIX 5: Create a local NumPy Generator.
    # # You aren't currently doing random operations in NumPy here, but if you
    # # ever add `rng.shuffle(X_train_np)`, this ensures it stays deterministic.
    # rng = np.random.default_rng(seed)

    objective = "reg:absoluteerror" if loss_type == "L1" else "reg:squarederror"
    print(f"--- Starting Pipeline | Objective: {objective} | Eval: {eval_group} | Seed: {seed} ---")

    # <-- FIX 6: Explicitly sort the dataframe before doing anything to prevent row-order shifting
    final_data = final_data.sort_values(by=features).reset_index(drop=True)

    train_df, val_df, test_df = split_data(
        final_data,
        features,
        split_col=split_config["col"],
        train_val=split_config["train"],
        val_val=split_config["val"],
        test_val=split_config["test"],
    )

    X_train_np, y_train_np = train_df[features].values, train_df[INHIBITION].values
    X_val_np, y_val_np = val_df[features].values, val_df[INHIBITION].values
    X_test_np, y_test_np = test_df[features].values, test_df[INHIBITION].values

    X_train_gpu, y_train_gpu = cp.array(X_train_np), cp.array(y_train_np)
    X_val_gpu, X_test_gpu = cp.array(X_val_np), cp.array(X_test_np)

    train_eval_idx = get_large_cohort_indices(train_df, eval_group)
    val_eval_idx = get_large_cohort_indices(val_df, eval_group)
    test_eval_idx = get_large_cohort_indices(test_df, eval_group)
    val_select_idx = get_large_cohort_indices(val_df, select_group)
    test_select_idx = get_large_cohort_indices(test_df, select_group)

    val_eval_groups_optuna = []
    for _, group in val_df.reset_index(drop=True).groupby(eval_group):
        if len(group) >= 20:
            val_eval_groups_optuna.append((group.index.values, y_val_np[group.index.values]))

    print("\n[1/3] Running Optuna Hyperparameter Tuning...")
    # <-- Pass seed here
    best_params = tune_hyperparameters(
        train_data=(X_train_np, y_train_np),
        val_data=(X_val_np, y_val_np),
        objective=objective,
        val_eval_groups_optuna=val_eval_groups_optuna,
        seed=seed,
    )
    print(f"Tuning Complete. Best Params: {best_params}")

    print("\n[2/3] Running Backward Feature Selection...")
    # <-- Pass seed here
    optimal_features, selection_history_df = run_backward_selection(
        train_data=(X_train_gpu, y_train_gpu, y_train_np),
        val_data=(X_val_gpu, y_val_np),
        test_data=(X_test_gpu, y_test_np),
        features=features,
        best_params=best_params,
        val_eval_idx=val_eval_idx,
        val_select_idx=val_select_idx,
        test_select_idx=test_select_idx,
        seed=seed,
    )
    print(f"\nOptimal Parsimonious Features Selected: {len(optimal_features)}")

    print("\n[3/3] Training Final Parsimonious Models...")

    X_train_val_gpu = cp.vstack((X_train_gpu, X_val_gpu))
    y_train_val_gpu = cp.concatenate((y_train_gpu, cp.array(y_val_np)))

    print("  -> Training Model 1: Train + Val Data")
    # <-- Pass seed here
    model_train_val = train_final_model(
        optimal_features, features, train_data=(X_train_val_gpu, y_train_val_gpu), best_params=best_params, seed=seed
    )

    X_all_gpu = cp.vstack((X_train_val_gpu, X_test_gpu))
    y_all_gpu = cp.concatenate((y_train_val_gpu, cp.array(y_test_np)))

    print("  -> Training Model 2: Entire Dataset (Train + Val + Test)")
    # <-- Pass seed here
    model_all = train_final_model(
        optimal_features, features, train_data=(X_all_gpu, y_all_gpu), best_params=best_params, seed=seed
    )

    # Evaluate & Save logic...
    feat_to_idx = {f: i for i, f in enumerate(features)}
    opt_idxs = [feat_to_idx[f] for f in optimal_features]
    X_train_opt = X_train_np[:, opt_idxs]
    X_val_opt = X_val_np[:, opt_idxs]

    # 2. Calculate Metrics for Model 1 (Train + Val)
    metrics_train_val = {
        "Train": evaluate_final_performance(model_train_val, X_train_opt, y_train_np, train_eval_idx, optimal_features),
        "Validation": evaluate_final_performance(model_train_val, X_val_opt, y_val_np, val_eval_idx, optimal_features),
    }

    # 3. Calculate Metrics for Model 2 (All Data)
    metrics_all = {
        "Train": evaluate_final_performance(model_all, X_train_opt, y_train_np, train_eval_idx, optimal_features),
        "Validation": evaluate_final_performance(model_all, X_val_opt, y_val_np, val_eval_idx, optimal_features),
    }

    # 4. Define Output Paths
    path_train_val = save_path.replace(".json", "_TrainVal.json")
    path_all = save_path.replace(".json", "_AllData.json")
    path_history = save_path.replace(".json", "_RFE_History.csv")

    # 5. Save XGBoost Models
    model_train_val.save_model(path_train_val)
    model_all.save_model(path_all)

    # 6. Save Metrics to JSON
    with open(path_train_val.replace(".json", "_metrics.json"), "w") as f:
        json.dump(metrics_train_val, f, indent=4)

    with open(path_all.replace(".json", "_metrics.json"), "w") as f:
        json.dump(metrics_all, f, indent=4)

    # 7. Save RFE Selection History to CSV
    selection_history_df.to_csv(path_history, index=False)

    print(f"\nPipeline Complete.")
    print(f"  - Model (Train+Val) saved to : {path_train_val}")
    print(f"  - Metrics saved to           : {path_train_val.replace('.json', '_metrics.json')}")
    print(f"  - Model (All Data) saved to  : {path_all}")
    print(f"  - RFE History saved to       : {path_history}")

    return {'model_train_val': model_train_val, 'model_all': model_all}, optimal_features, selection_history_df