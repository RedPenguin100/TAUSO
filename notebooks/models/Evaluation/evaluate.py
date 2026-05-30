"""
Evaluate the TAUSO model against the competitors on the test split (Part 1 of the pipeline).

Writes tidy artifacts under results/ that the graphs notebooks
(notebooks/graphs/models/) consume directly. Re-run whenever the model, data, or competitor
set changes and every downstream graph regenerates unchanged.

Scorers (each a single prediction column scored against Inhibition(%)):
  TAUSO               our model's prediction on the held-out test rows
  oligo_ai_score      OligoAI (the key, calibrated competitor)
  PFRED_PLS, PFRED_SVM | OW_Overall, OW_Tm, OW_Intra_Oligo, OW_Duplex | sfold_accessibility |
  miranda_score, miranda_energy

Cohorts:
  patent            one patent inhibition table = one experiment (custom_id, the per-experiment
                    ASO-selection unit; nested inside patent_id)
  gene x cell_line  biological cohort (canonical gene x cell line)

Breakdowns (metrics recomputed within each slice, with per-slice n and inhibition range):
  chemistry   cEt vs 2'MOE gapmers
  cell_line   per cell line

Honest-comparison rules (see the compare-to-oligoai skill): competitors are NOT model inputs;
each scorer is compared only where present; competitor sign is oriented on the TRAIN split;
RMSE/MAE are reported only for the two models that predict calibrated inhibition (TAUSO,
OligoAI) -- the others are rankers/biophysical scores and absolute error is not their goal.
OligoAI sanity floor: median per-cohort Spearman >=~0.40, else an index_oligo alignment bug.

Usage:
    python -u -m notebooks.models.Evaluation.evaluate
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
CALIBRATED = {MODEL_NAME, "oligo_ai_score"}        # only these aim at absolute %inhibition
COMPETITORS = [
    "oligo_ai_score",
    "PFRED_PLS", "PFRED_SVM",
    "OW_Overall", "OW_Tm", "OW_Intra_Oligo", "OW_Duplex",
    "sfold_accessibility",
    "miranda_score", "miranda_energy",
]
COHORT = "custom_id"                                # the per-experiment unit
GROUPINGS = {"patent": [COHORT], "gene x cell_line": [CANONICAL_GENE, CELL_LINE]}
DEFAULT_MODEL = Path("notebooks/models/Model_Oligo_TrainVal.json")
DEFAULT_OUT = Path(__file__).parent / "results"
MIN_COHORT = 5
OLIGOAI_FLOOR = 0.30


def _chem_label(mod):
    if not isinstance(mod, str):
        return "other"
    return "cEt" if "cEt" in mod else ("2'MOE" if "MOE" in mod else "other")


# --------------------------------------------------------------------------- metrics

def _group_metrics(t, p):
    """Per-cohort metrics for one group's (truth, score) arrays, or None if unscoreable."""
    finite = np.isfinite(t) & np.isfinite(p)
    t, p = t[finite], p[finite]
    n = len(t)
    if n < MIN_COHORT or np.ptp(t) == 0 or np.ptp(p) == 0:
        return None
    corr, _ = spearmanr(t, p)
    if np.isnan(corr):
        return None
    out = {"n": n, "spearman": float(corr)}
    rel, scr = np.clip(t, 0, None).reshape(1, -1), p.reshape(1, -1)
    for k in (5, 10):
        kk = min(k, n)
        top_pred = np.argpartition(p, -kk)[-kk:]
        out[f"top{k}_inhib"] = float(np.mean(t[top_pred]))
        if n > k:
            out[f"ndcg{k}"] = float(ndcg_score(rel, scr, k=k))
            out[f"p@{k}"] = np.intersect1d(top_pred, np.argpartition(t, -kk)[-kk:]).size / kk
    return out


def _aggregate(per, scorer, **extra):
    agg = {"scorer": scorer, "n_cohorts": len(per), **extra}
    for m in ("spearman", "top5_inhib", "top10_inhib", "ndcg5", "ndcg10", "p@5", "p@10"):
        vals = [d[m] for d in per if m in d]
        if vals:
            agg[f"{m}_median"], agg[f"{m}_mean"] = float(np.median(vals)), float(np.mean(vals))
    return agg


def _per_cohort(df, score, cohort_keys):
    """List of per-cohort metric dicts for one score array over the given integer cohort keys."""
    t = df[INHIBITION].to_numpy("float64")
    per = []
    for g in np.unique(cohort_keys):
        res = _group_metrics(t[cohort_keys == g], score[cohort_keys == g])
        if res is not None:
            per.append(res)
    return per


def evaluate(preds_df, scorers):
    """Returns (summary_df, per_group_df) over the headline cohort groupings."""
    summary_rows, per_group_rows = [], []
    for scorer in scorers:
        score = preds_df[scorer].to_numpy("float64")
        for gname, gcols in GROUPINGS.items():
            keys = preds_df.groupby(gcols).ngroup().to_numpy()
            per = _per_cohort(preds_df, score, keys)
            if not per:
                continue
            summary_rows.append(_aggregate(per, scorer, grouping=gname))
            per_group_rows += [{"scorer": scorer, "grouping": gname,
                                "spearman": d["spearman"], "n": d["n"]} for d in per]
        # global calibration (RMSE/MAE only for the calibrated models)
        t, p = preds_df[INHIBITION].to_numpy("float64"), score
        ok = np.isfinite(t) & np.isfinite(p)
        if ok.sum() >= 2:
            row = {"scorer": scorer, "grouping": "GLOBAL", "n_cohorts": 1,
                   "global_spearman": float(spearmanr(t[ok], p[ok])[0])}
            if scorer in CALIBRATED:
                row["MAE"] = float(mean_absolute_error(t[ok], p[ok]))
                row["RMSE"] = float(np.sqrt(mean_squared_error(t[ok], p[ok])))
            summary_rows.append(row)
    return pd.DataFrame(summary_rows), pd.DataFrame(per_group_rows)


def breakdown(preds_df, scorers, slice_col):
    """Metrics per slice of `slice_col`, cohorts = custom_id within the slice, plus per-slice
    descriptive stats (n samples, inhibition range)."""
    rows = []
    for sval, sub in preds_df.groupby(slice_col):
        keys = sub.groupby(COHORT).ngroup().to_numpy()
        inh = sub[INHIBITION]
        desc = {"breakdown": slice_col, "slice": sval, "n_samples": len(sub),
                "n_genes": int(sub[CANONICAL_GENE].nunique()), "n_patents": int(sub[COHORT].nunique()),
                "inhib_min": float(inh.min()), "inhib_max": float(inh.max()),
                "inhib_mean": float(inh.mean()), "base_hit70": float((inh >= 70).mean())}
        for scorer in scorers:
            per = _per_cohort(sub, sub[scorer].to_numpy("float64"), keys)
            if per:
                rows.append(_aggregate(per, scorer, **desc))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- driver

def build_predictions(model_path):
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
    present = [c for c in COMPETITORS if c in final_data.columns]
    missing = [c for c in COMPETITORS if c not in final_data.columns]
    if missing:
        logger.warning("competitors missing (skipped): %s", missing)
    scorers = [MODEL_NAME] + present

    test_df[MODEL_NAME] = bst.predict(xgb.DMatrix(test_df[feats].to_numpy("float32"), feature_names=feats))
    test_df["chemistry"] = test_df["Modification"].map(_chem_label)

    # orient competitor signs on the TRAIN split (higher == more inhibition)
    y = train_df[INHIBITION].to_numpy("float64")
    tk = train_df.groupby(COHORT).ngroup().to_numpy()
    for c in present:
        per = _per_cohort(train_df.assign(_s=train_df[c]), train_df[c].to_numpy("float64"), tk)
        med = np.median([d["spearman"] for d in per]) if per else 0.0
        if med < 0:
            test_df[c] = -test_df[c].to_numpy("float64")
            logger.info("Oriented '%s' with sign -1 (train median Spearman %.3f)", c, med)
    _ = y

    keep = ["index_oligo", INHIBITION, "patent_id", "custom_id", "chemistry",
            CANONICAL_GENE, CELL_LINE] + scorers
    keep = [c for c in keep if c in test_df.columns]
    logger.info("test rows: %d | scorers: %d", len(test_df), len(scorers))
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
    chem_df = breakdown(preds_df, scorers, "chemistry")
    cell_df = breakdown(preds_df, scorers, CELL_LINE)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    preds_df.to_parquet(args.out_dir / "test_predictions.parquet", index=False)
    summary_df.to_csv(args.out_dir / "metrics_summary.csv", index=False)
    per_group_df.to_csv(args.out_dir / "per_group_spearman.csv", index=False)
    pd.concat([chem_df, cell_df], ignore_index=True).to_csv(args.out_dir / "breakdown_metrics.csv", index=False)

    sp = summary_df[summary_df.grouping == "patent"].set_index("scorer")["spearman_median"].round(3)
    logger.info("Median per-patent-table Spearman:\n%s", sp.sort_values(ascending=False).to_string())
    oa = sp.get("oligo_ai_score", float("nan"))
    if np.isfinite(oa) and oa < OLIGOAI_FLOOR:
        logger.warning("OligoAI median Spearman %.3f < %.2f floor — check index_oligo alignment!", oa, OLIGOAI_FLOOR)
    logger.info("Saved -> %s/{test_predictions.parquet,metrics_summary.csv,per_group_spearman.csv,breakdown_metrics.csv}", args.out_dir)


if __name__ == "__main__":
    main()
