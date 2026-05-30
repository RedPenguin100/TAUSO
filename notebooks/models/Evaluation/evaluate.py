"""
Evaluate the TAUSO model against the competitors on the test split (Part 1 of the pipeline).

Produces three tidy artifacts under results/ that the graphs notebook
(notebooks/graphs/models/model_comparison.ipynb) consumes directly. Re-run this whenever the
model, data, or competitor set changes and every downstream graph regenerates unchanged.

Scorers compared (each a single prediction column scored against Inhibition(%)):
  TAUSO               our model's prediction on the held-out test rows
  oligo_ai_score      OligoAI (the key competitor)
  PFRED_PLS, PFRED_SVM
  OW_Overall, OW_Tm, OW_Intra_Oligo, OW_Duplex   (OligoWalk)
  sfold_accessibility (sfold)
  miranda_score, miranda_energy

Honest-comparison rules (see the compare-to-oligoai skill):
  * competitors are NOT in the model's input features (model trained load_competition=False);
  * each scorer is compared only on test rows where that scorer is present;
  * competitor directionality is oriented on the TRAIN split (flip sign if train per-cohort
    median Spearman is negative) so accessibility/energy scores are not unfairly penalised;
  * OligoAI sanity floor: its median per-cohort Spearman should be >=~0.40 — below that is an
    index_oligo alignment bug, not a real result.

Groupings (every metric is computed per group, then aggregated):
  patent_id            the patent each oligo came from
  gene x cell_line     biological cohort (canonical gene x cell line)
  custom_id            within-experiment cohort (the ASO-selection setting)

Usage:
    python -u -m notebooks.models.Evaluation.evaluate
    python -u -m notebooks.models.Evaluation.evaluate --model notebooks/models/Model_Oligo_TrainVal.json
"""
import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, ndcg_score

from notebooks.models.utility import load_and_validate_final_data
from tauso.data.consts import CANONICAL_GENE, CELL_LINE, INHIBITION

logger = logging.getLogger(__name__)

MODEL_NAME = "TAUSO"
COMPETITORS = [
    "oligo_ai_score",
    "PFRED_PLS", "PFRED_SVM",
    "OW_Overall", "OW_Tm", "OW_Intra_Oligo", "OW_Duplex",
    "sfold_accessibility",
    "miranda_score", "miranda_energy",
]
GROUPINGS = {
    "patent_id":      ["patent_id"],
    "gene x cell_line": [CANONICAL_GENE, CELL_LINE],
    "custom_id":      ["custom_id"],
}
DEFAULT_MODEL = Path("notebooks/models/Model_Oligo_TrainVal.json")
DEFAULT_OUT = Path(__file__).parent / "results"
MIN_COHORT = 5          # min rows in a group to score it (need a rank to correlate)
OLIGOAI_FLOOR = 0.30    # warn if OligoAI median per-cohort Spearman drops below this


# --------------------------------------------------------------------------- metrics

def _group_metrics(t, p):
    """All per-cohort metrics for one group's (truth, score) arrays, or None if unscoreable."""
    finite = np.isfinite(t) & np.isfinite(p)
    t, p = t[finite], p[finite]
    n = len(t)
    if n < MIN_COHORT or np.ptp(t) == 0 or np.ptp(p) == 0:
        return None
    corr, _ = spearmanr(t, p)
    if np.isnan(corr):
        return None
    out = {"n": n, "spearman": float(corr)}
    rel = np.clip(t, 0, None).reshape(1, -1)
    scr = p.reshape(1, -1)
    for k in (5, 10):
        kk = min(k, n)
        top_pred = np.argpartition(p, -kk)[-kk:]
        out[f"top{k}_inhib"] = float(np.mean(t[top_pred]))          # mean true inhib of model's best k
        if n > k:
            out[f"ndcg{k}"] = float(ndcg_score(rel, scr, k=k))
            top_true = np.argpartition(t, -kk)[-kk:]
            out[f"p@{k}"] = np.intersect1d(top_pred, top_true).size / kk
    return out


def _orient_signs(train_df, scorers):
    """Sign per scorer that makes higher == more inhibition, decided on the TRAIN split."""
    signs = {}
    y = train_df[INHIBITION].to_numpy("float64")
    keys = train_df.groupby("custom_id").ngroup().to_numpy()
    for s in scorers:
        if s == MODEL_NAME:
            signs[s] = 1.0
            continue
        vals = train_df[s].to_numpy("float64")
        corrs = []
        for g in np.unique(keys):
            m = keys == g
            t, p = y[m], vals[m]
            ok = np.isfinite(t) & np.isfinite(p)
            if ok.sum() >= MIN_COHORT and np.ptp(t[ok]) and np.ptp(p[ok]):
                c, _ = spearmanr(t[ok], p[ok])
                if not np.isnan(c):
                    corrs.append(c)
        med = np.median(corrs) if corrs else 0.0
        signs[s] = -1.0 if med < 0 else 1.0
        if signs[s] < 0:
            logger.info("Oriented '%s' with sign -1 (train median Spearman %.3f)", s, med)
    return signs


def evaluate(preds_df, scorers):
    """preds_df: one row per test oligo with INHIBITION + grouping cols + one column per scorer.

    Returns (summary_df, per_group_df):
      summary_df   : scorer x grouping x metric -> aggregate (median + mean over groups, global)
      per_group_df : scorer x grouping x group_id -> per-cohort spearman/n (for box plots)
    """
    summary_rows, per_group_rows = [], []
    for scorer in scorers:
        score_all = preds_df[scorer].to_numpy("float64")
        for gname, gcols in GROUPINGS.items():
            keys = preds_df.groupby(gcols).ngroup().to_numpy()
            per = []
            for g in np.unique(keys):
                m = keys == g
                res = _group_metrics(preds_df[INHIBITION].to_numpy("float64")[m], score_all[m])
                if res is None:
                    continue
                per.append(res)
                per_group_rows.append({"scorer": scorer, "grouping": gname,
                                       "spearman": res["spearman"], "n": res["n"]})
            if not per:
                continue
            agg = {"scorer": scorer, "grouping": gname, "n_cohorts": len(per)}
            for metric in ("spearman", "top5_inhib", "top10_inhib", "ndcg5", "ndcg10", "p@5", "p@10"):
                vals = [d[metric] for d in per if metric in d]
                if vals:
                    agg[f"{metric}_median"] = float(np.median(vals))
                    agg[f"{metric}_mean"] = float(np.mean(vals))
            summary_rows.append(agg)

    # global (un-grouped) calibration metrics
    for scorer in scorers:
        t = preds_df[INHIBITION].to_numpy("float64")
        p = preds_df[scorer].to_numpy("float64")
        ok = np.isfinite(t) & np.isfinite(p)
        if ok.sum() < 2:
            continue
        gc, _ = spearmanr(t[ok], p[ok])
        summary_rows.append({"scorer": scorer, "grouping": "GLOBAL", "n_cohorts": 1,
                             "global_spearman": float(gc),
                             "MAE": float(mean_absolute_error(t[ok], p[ok])),
                             "RMSE": float(np.sqrt(mean_squared_error(t[ok], p[ok])))})
    return pd.DataFrame(summary_rows), pd.DataFrame(per_group_rows)


# --------------------------------------------------------------------------- driver

def build_predictions(model_path):
    """Load data + model, return the test-set predictions table (one column per scorer)."""
    logger.info("Loading data (with competition baselines)...")
    final_data, _ = load_and_validate_final_data(version="oligo", load_competition=True)
    bst = xgb.Booster()
    bst.load_model(str(model_path))
    feats = bst.feature_names
    leaked = [c for c in COMPETITORS if c in feats]
    if leaked:
        raise SystemExit(f"Model uses competitor columns as features: {leaked}. Retrain with load_competition=False.")

    train_df = final_data[final_data["split"] == "train"]
    test_df = final_data[final_data["split"] == "test"].copy()
    logger.info("test rows: %d | model features: %d", len(test_df), len(feats))

    present = [c for c in COMPETITORS if c in final_data.columns]
    missing = [c for c in COMPETITORS if c not in final_data.columns]
    if missing:
        logger.warning("competitors missing (skipped): %s", missing)
    scorers = [MODEL_NAME] + present

    test_df[MODEL_NAME] = bst.predict(xgb.DMatrix(test_df[feats].to_numpy("float32"), feature_names=feats))
    signs = _orient_signs(train_df, present)
    for c in present:
        test_df[c] = test_df[c].to_numpy("float64") * signs[c]

    keep = [INHIBITION, "patent_id", CANONICAL_GENE, CELL_LINE, "custom_id"] + scorers
    if "index_oligo" in test_df.columns:
        keep = ["index_oligo"] + keep
    return test_df[keep].reset_index(drop=True), scorers


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s %(levelname)s | %(message)s",
                        stream=sys.stdout, force=True)

    preds_df, scorers = build_predictions(args.model)
    summary_df, per_group_df = evaluate(preds_df, scorers)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    preds_df.to_parquet(args.out_dir / "test_predictions.parquet", index=False)
    summary_df.to_csv(args.out_dir / "metrics_summary.csv", index=False)
    per_group_df.to_csv(args.out_dir / "per_group_spearman.csv", index=False)

    # headline + sanity check
    sp = (summary_df[summary_df.grouping == "custom_id"]
          .set_index("scorer")["spearman_median"].round(3))
    logger.info("Median per-cohort (custom_id) Spearman:\n%s", sp.sort_values(ascending=False).to_string())
    oa = sp.get("oligo_ai_score", float("nan"))
    if np.isfinite(oa) and oa < OLIGOAI_FLOOR:
        logger.warning("OligoAI median Spearman %.3f < %.2f floor — check index_oligo alignment!", oa, OLIGOAI_FLOOR)
    logger.info("Saved -> %s/{test_predictions.parquet,metrics_summary.csv,per_group_spearman.csv}", args.out_dir)


if __name__ == "__main__":
    main()
