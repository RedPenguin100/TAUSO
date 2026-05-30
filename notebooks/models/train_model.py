"""
Train the single TAUSO ASO-efficacy model: a plain XGBoost regressor on train+val.

Part 1 of the article pipeline. One model, fixed regularized params, no grid, no feature
selection. Features come strictly from the branch's wide feature cache (load_competition=False
keeps every competitor score out of the inputs, so the head-to-head in Evaluation is honest).

The model is saved with its feature_names baked in so the evaluator reads them directly:
  <out>/Model_Oligo_TrainVal.json   trained on train+val (test held out — use this for evaluation)

Usage:
    python -u -m notebooks.models.train_model
    python -u -m notebooks.models.train_model --out notebooks/models/Model_Oligo_TrainVal.json
"""
import argparse
import logging
import sys
from pathlib import Path

import xgboost as xgb

from notebooks.models.utility import load_and_validate_final_data
from tauso.data.consts import INHIBITION

logger = logging.getLogger(__name__)

DEFAULT_OUT = Path(__file__).parent / "Model_Oligo_TrainVal.json"

# Fixed, regularized params. Plain L2 regression on CPU.
PARAMS = {
    "objective": "reg:squarederror",
    "tree_method": "hist", "device": "cpu",
    "max_depth": 6, "learning_rate": 0.05,
    "subsample": 0.8, "colsample_bytree": 0.6, "min_child_weight": 20,
    "reg_lambda": 5.0, "reg_alpha": 1.0, "seed": 1,
}
ROUNDS = 700


def train(final_data, features):
    """Train one XGBoost regressor on the train+val rows; return the fitted booster."""
    trv = final_data[final_data["split"].isin(["train", "val"])]
    logger.info("Training on train+val: %d rows, %d features", len(trv), len(features))
    X = trv[features].to_numpy("float32")
    y = trv[INHIBITION].to_numpy("float64")
    dmat = xgb.QuantileDMatrix(X, label=y, feature_names=features)
    return xgb.train(PARAMS, dmat, num_boost_round=ROUNDS)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="Where to save the model JSON")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s %(levelname)s | %(message)s",
                        stream=sys.stdout, force=True)

    logger.info("Loading data (cache features only, no competition scores)...")
    final_data, features = load_and_validate_final_data(version="oligo", load_competition=False)
    logger.info("Loaded %d rows, %d candidate features", len(final_data), len(features))

    booster = train(final_data, features)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    booster.save_model(str(args.out))
    logger.info("Saved model -> %s", args.out)
    logger.info("Evaluate with: python -m notebooks.models.Evaluation.evaluate --model %s", args.out)


if __name__ == "__main__":
    main()
