"""Step 5 -- the final model. Hardcodes the chosen configuration: clean_exp objective, the top-K RMSE-ranked
features (K below, from 4_feature_selection.py's feature_ranking.txt), and the regularized params.

Two data modes:
  --data trainval : train on train+val, evaluate ONCE on the held-out test -- the reported performance.
  --data all      : train on ALL data (train+val+test) for deployment / scoring new candidates (e.g. DIAPH3);
                    no test evaluation (test is in the training set).

Saves the trained model (+ the feature list, + test predictions in trainval mode).

Run:  python notebooks/models/5_final_model.py [--data trainval|all]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # repo root, for the notebooks.* imports
from notebooks.models import common
from notebooks.models.evaluate import evaluate
from tauso.data.consts import INHIBITION_PERCENT

VARIANT = "clean_exp"
K = 450
PARAMS = dict(tree_method="hist", device="cuda", max_depth=6, learning_rate=0.05,
              subsample=0.8, colsample_bytree=0.6, min_child_weight=20, reg_lambda=5.0, reg_alpha=1.0)
ROUNDS = 700
SEEDS = (1, 2, 3)
METR = ["exp_med", "exp_mean", "gxc_med", "gxc_mean", "top5", "p5", "p10", "globP"]
OLIGOAI_CSV = "/home/michael/career/tauso_article/.tauso_data/features/oligo/_patches/oligo_ai_score.csv"


def chosen_cols(features):
    ranking = (common.RESULTS_DIR / f"feature_ranking_{VARIANT}.txt").read_text().splitlines()
    keep = set(ranking[:K])
    return [f for f in features if f in keep]                      # canonical order


def suite(pred, df):
    return evaluate(pred, df[INHIBITION_PERCENT].to_numpy(np.float64),
                    df["custom_id"].to_numpy(), df["cohort_id"].to_numpy())


def fmt(v, m):
    return f"{v:.1f}" if m == "top5" else f"{v:.3f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", choices=["trainval", "all"], default="trainval")
    args = ap.parse_args()

    df, features = common.load_dataset()
    cols = chosen_cols(features)
    trv, test = common.split(df)
    common.write_checked(common.RESULTS_DIR / f"final_features_K{K}.txt", "\n".join(cols) + "\n")

    if args.data == "all":
        model = common.train(df, cols, VARIANT, PARAMS, ROUNDS, seed=1)
        out = common.RESULTS_DIR / f"final_model_{VARIANT}_K{K}_all.json"
        model.save_model(str(out))
        print(f"trained {VARIANT} @ K={K} on ALL {len(df)} rows, {len(cols)} features -> {out}  (deployment; no test eval)")
        return

    # trainval: train on train+val, evaluate ONCE on the held-out test (3 seeds)
    runs = [suite(common.predict(common.train(trv, cols, VARIANT, PARAMS, ROUNDS, s), test, cols), test) for s in SEEDS]
    agg = {m: (np.mean([r[m] for r in runs]), np.std([r[m] for r in runs])) for m in METR}

    te = test.merge(pd.read_csv(OLIGOAI_CSV), on="index_oligo", how="left")
    cov = te["oligo_ai_score"].notna().mean()
    te_oai = te[te["oligo_ai_score"].notna()]
    oai = suite(te_oai["oligo_ai_score"].to_numpy(np.float64), te_oai)

    print(f"\nFINAL MODEL  {VARIANT} @ K={K}  |  train+val ({len(trv)}) -> held-out test ({len(test)})  |  {len(SEEDS)} seeds\n")
    print("           " + " ".join(f"{m:>13}" for m in METR))
    print("TAUSO      " + " ".join(f"{fmt(agg[m][0], m)}±{fmt(agg[m][1], m)}".rjust(13) for m in METR))
    print("OligoAI    " + " ".join(fmt(oai[m], m).rjust(13) for m in METR) + f"   (test cov {cov:.0%})")

    model = common.train(trv, cols, VARIANT, PARAMS, ROUNDS, seed=1)
    model.save_model(str(common.RESULTS_DIR / f"final_model_{VARIANT}_K{K}_trainval.json"))
    np.savez(common.RESULTS_DIR / f"final_test_pred_K{K}.npz",
             index_oligo=test["index_oligo"].to_numpy(), pred=common.predict(model, test, cols))
    print(f"\nsaved -> final_model_{VARIANT}_K{K}_trainval.json, final_test_pred_K{K}.npz, final_features_K{K}.txt")


if __name__ == "__main__":
    main()
