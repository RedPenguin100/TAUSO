"""
Train simple XGBoost ASO-efficacy models — fixed params, no feature selection.

Trains fixed-param models over a small, interpretable grid; the evaluation
(Evaluation/evaluate_all.py) reports which config wins on which metric / chemistry / cell line.

Grid axes:
  loss  : L2 (reg:squarederror) | L1 (reg:absoluteerror) | ndcg (rank:ndcg)
  rbp   : norbp (drop RBP_*) | withrbp (keep them)
  sense : both | raw (drop normalized n_sense_*) | n (drop raw sense_*, keep n_sense_*)
  split : (NDCG only) cohort=gene×cell-line | custom_id  -- the query-group for ranking.
          Regression is loss-on-rows (no grouping), so it is split-agnostic and trained once;
          evaluation reports it on both groupings anyway.

Each model is saved with feature_names baked in so the evaluator reads them directly:
  <out-dir>/<config>/Model_Oligo_<config>_TrainVal.json   (trained on train+val)
  <out-dir>/<config>/Model_Oligo_<config>_AllData.json    (train+val+test; deployment, test-leaked)

Usage:
    python -u -m notebooks.models.train_model                      # full grid
    python -u -m notebooks.models.train_model --configs L2_norbp_both,L1_norbp_both
    python -u -m notebooks.models.train_model --list
"""
import argparse
import itertools
import logging
import sys
from pathlib import Path

import numpy as np
import xgboost as xgb

from notebooks.models.utility import load_and_validate_final_data
from tauso.data.consts import CANONICAL_GENE, CELL_LINE, INHIBITION

logger = logging.getLogger(__name__)

DEFAULT_OUTPUT_DIR = Path(__file__).parent / "SimpleModels"

# Fixed, sensible, regularized params (the config that beat the tuned pipeline).
SANE = {
    "tree_method": "hist", "device": "cuda",
    "max_depth": 6, "learning_rate": 0.05,
    "subsample": 0.8, "colsample_bytree": 0.6, "min_child_weight": 20,
    "reg_lambda": 5.0, "reg_alpha": 1.0, "seed": 1,
}
ROUNDS = 700
OBJECTIVES = {"L2": "reg:squarederror", "L1": "reg:absoluteerror", "ndcg": "rank:ndcg"}
RAW_SENSE_DROP = r"^sense_(3utr|5utr|exon|intron|utr)$"  # raw versions that have an n_ counterpart


# ---------------------------------------------------------------------------
# Feature sets
# ---------------------------------------------------------------------------

def _feature_set(all_features, rbp, sense):
    import re
    feats = list(all_features)
    if rbp == "norbp":
        feats = [f for f in feats if not f.startswith("RBP_")]
    if sense == "raw":                                   # keep raw sense_*, drop normalized n_sense_*
        feats = [f for f in feats if not f.startswith("n_")]
    elif sense == "n":                                   # keep n_sense_*, drop raw sense_3utr/.../utr
        feats = [f for f in feats if not re.search(RAW_SENSE_DROP, f)]
    return feats


# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------

def build_grid():
    """Return {config_name: dict(loss, rbp, sense, split)}."""
    grid = {}
    rbps, senses = ["norbp", "withrbp"], ["both", "raw", "n"]
    for loss in ("L2", "L1"):                            # regression: split-agnostic
        for rbp, sense in itertools.product(rbps, senses):
            grid[f"{loss}_{rbp}_{sense}"] = dict(loss=loss, rbp=rbp, sense=sense, split=None)
    for split in ("cohort", "custom_id"):                # ndcg: query group depends on split
        for rbp, sense in itertools.product(rbps, senses):
            grid[f"ndcg_{rbp}_{sense}_{split}"] = dict(loss="ndcg", rbp=rbp, sense=sense, split=split)
    return grid


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------

def _ndcg_groups(df, split):
    """Contiguous-group sort order + group sizes for a rank:ndcg query grouping."""
    key = df.groupby([CANONICAL_GENE, CELL_LINE]).ngroup() if split == "cohort" else df["custom_id"].astype("category").cat.codes
    g = key.to_numpy()
    order = np.argsort(g, kind="stable")
    sizes = np.bincount(g[order])
    return order, sizes[sizes > 0]


def _train(df, feats, cfg):
    y = df[INHIBITION].to_numpy(np.float64)
    X = df[feats].to_numpy(np.float32)
    params = {**SANE, "objective": OBJECTIVES[cfg["loss"]]}
    if cfg["loss"] == "ndcg":
        order, sizes = _ndcg_groups(df, cfg["split"])
        grades = np.rint(np.clip(y, 0, 100) / 100.0 * 15.0).astype(np.int32)[order]  # 0..15 relevance
        d = xgb.DMatrix(X[order], label=grades, feature_names=feats)
        d.set_group(sizes)
    else:
        d = xgb.QuantileDMatrix(X, label=y, feature_names=feats)
    return xgb.train(params, d, num_boost_round=ROUNDS)


def train_config(name, cfg, final_data, all_features, out_dir):
    feats = _feature_set(all_features, cfg["rbp"], cfg["sense"])
    trv = final_data[final_data["split"].isin(["train", "val"])]
    alld = final_data
    logger.info("[%s] loss=%s rbp=%s sense=%s split=%s | %d features",
                name, cfg["loss"], cfg["rbp"], cfg["sense"], cfg["split"], len(feats))

    cdir = out_dir / name
    cdir.mkdir(parents=True, exist_ok=True)
    _train(trv,  feats, cfg).save_model(str(cdir / f"Model_Oligo_{name}_TrainVal.json"))
    _train(alld, feats, cfg).save_model(str(cdir / f"Model_Oligo_{name}_AllData.json"))
    logger.info("[%s] saved -> %s", name, cdir)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--configs", default=None, help="Comma-separated config names (default: full grid)")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--list", action="store_true", help="List grid configs and exit")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s %(levelname)s | %(message)s",
                        stream=sys.stdout, force=True)

    grid = build_grid()
    if args.list:
        print(f"\n{len(grid)} configs:")
        for n, c in grid.items():
            print(f"  {n}")
        return

    selected = args.configs.split(",") if args.configs else list(grid)
    unknown = [c for c in selected if c not in grid]
    if unknown:
        logger.error("Unknown configs: %s", unknown); sys.exit(1)

    logger.info("Loading data...")
    final_data, all_features = load_and_validate_final_data(version="oligo")
    logger.info("Training %d configs -> %s", len(selected), args.output)
    for name in selected:
        train_config(name, grid[name], final_data, all_features, args.output)
    logger.info("Done. Evaluate with: python -m notebooks.models.Evaluation.evaluate_all --models-dir %s", args.output)


if __name__ == "__main__":
    main()
