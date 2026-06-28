"""Shared helpers for the numbered model-selection scripts (steps 1-5):
dataset loading, the frozen train/test split, frozen gene-grouped CV folds, evaluation, and reporting.
Keeps each script thin and consistent, and the split/CV definition in exactly one place.
"""
from pathlib import Path

import numpy as np
import xgboost as xgb
from sklearn.model_selection import GroupKFold

from notebooks.models.utility import load_and_validate_final_data
from notebooks.models.evaluate import evaluate
from notebooks.preprocessing import assign_cohort
from tauso.data.consts import CANONICAL_GENE_NAME, INHIBITION_PERCENT

METRICS = ("exp_med", "exp_mean", "gxc_med", "gxc_mean", "top5", "p5", "p10", "globP")
RESULTS_DIR = Path(__file__).parent / "results"


def load_dataset():
    """Canonical features + labels with the cohort_id column assigned. Returns (df, feature_names)."""
    df, features = load_and_validate_final_data(version="oligo")
    return assign_cohort(df), features


def split(df):
    """The frozen OligoAI split -> (train+val dataframe, test dataframe)."""
    return df[df["split"].isin(["train", "val"])], df[df["split"] == "test"]


def strict_gapmer_mask(df):
    """Boolean mask of strict gapmers: 5-10-5 2'-MOE or 3-10-3 cEt sugar pattern (STRICT_GAPMER_PATTERNS)
    AND a full phosphorothioate backbone (no PO linkages). Stricter than the repo's sugar-only filter."""
    from notebooks.preprocessing import STRICT_GAPMER_PATTERNS
    from tauso.data.consts import CHEMICAL_PATTERN
    return df[CHEMICAL_PATTERN].isin(STRICT_GAPMER_PATTERNS) & (df["mod_ps_po_percentage"] == 0)


def gene_cv_folds(trv, n_splits=5):
    """Frozen gene-grouped CV over train+val: yields (train_df, val_df) with no gene shared between them.
    Deterministic (GroupKFold uses no RNG), so the folds are identical on every run."""
    for tr, va in GroupKFold(n_splits).split(trv, groups=trv[CANONICAL_GENE_NAME]):
        yield trv.iloc[tr], trv.iloc[va]


def predict(model, df, features):
    return model.predict(xgb.DMatrix(df[features].to_numpy(np.float64), feature_names=features))


def metrics_on(model, df, features):
    """Full metric suite for `model` on `df`, ranked against actual inhibition."""
    return evaluate(predict(model, df, features), df[INHIBITION_PERCENT].to_numpy(np.float64),
                    df["custom_id"].to_numpy(), df["cohort_id"].to_numpy())


# variant -> (objective, query-group column, demean-by column, extra params)
VARIANTS = {
    "reg":           ("reg:squarederror", None,        None,        {}),
    "clean_exp":     ("reg:squarederror", None,        "custom_id", {}),
    "clean_gene":    ("reg:squarederror", None,        "cohort_id", {}),
    "pair_exp":      ("rank:pairwise",    "custom_id", None,        {}),
    "pair_gene":     ("rank:pairwise",    "cohort_id", None,        {}),
    "ndcg_exp":      ("rank:ndcg",        "custom_id", None,        {}),                       # exponential gain (default)
    "ndcg_gene":     ("rank:ndcg",        "cohort_id", None,        {}),
    "ndcg_exp_lin":  ("rank:ndcg",        "custom_id", None,        {"ndcg_exp_gain": False}),  # linear gain
    "ndcg_gene_lin": ("rank:ndcg",        "cohort_id", None,        {"ndcg_exp_gain": False}),
}


def _target(df, demean_by):
    y = df[INHIBITION_PERCENT].to_numpy(np.float64)
    if demean_by is not None:                                      # within-group deviation (the clean_* models)
        y = y - df.groupby(demean_by)[INHIBITION_PERCENT].transform("mean").to_numpy()
    return y


def train(df, features, variant, params, rounds, seed=1):
    """Train one variant (objective + target grouping from VARIANTS) with the given params on `df`."""
    objective, qid, demean_by, extra = VARIANTS[variant]
    config = {**params, "objective": objective, "seed": seed, **extra}
    if qid is not None:                                            # learning-to-rank: contiguous query groups
        df = df.sort_values(qid, kind="stable")
        label = df[INHIBITION_PERCENT].to_numpy(np.float64)
        if objective == "rank:ndcg":                               # NDCG ranks by integer relevance grades
            label = np.rint(np.clip(label, 0, 100) / 100.0 * 15.0)  # inhibition 0..100% -> grade 0..15
        dtrain = xgb.DMatrix(df[features].to_numpy(np.float64), label=label, feature_names=features)
        dtrain.set_info(group=df.groupby(qid, sort=False).size().to_numpy())
    else:
        dtrain = xgb.DMatrix(df[features].to_numpy(np.float64), label=_target(df, demean_by), feature_names=features)
    return xgb.train(config, dtrain, rounds)


def mean_std(scores_list):
    """{metric: (mean, std)} across a list of score dicts (CV folds or seeds)."""
    return {m: (float(np.mean([s[m] for s in scores_list])),
                float(np.std([s[m] for s in scores_list]))) for m in METRICS}


def metric_table(title, rows):
    """Aligned text table as a list of lines. rows = {label: {metric: value or (mean, std)}}."""
    out = [title, " " * 13 + " ".join(f"{m:>13}" for m in METRICS)]
    for label, scores in rows.items():
        cells = []
        for m in METRICS:
            v = scores[m]
            cell = f"{v[0]:.3f}±{v[1]:.3f}" if isinstance(v, tuple) else f"{v:.3f}"
            cells.append(f"{cell:>13}")
        out.append(f"{label:13s}" + " ".join(cells))
    return out


def write_checked(path, text):
    """Write `text` to `path`, warning first if a *different* version is already on disk.
    A cheap half-test of reproducibility: a deterministic re-run reproduces the file exactly, so a
    warning means this run's output changed -- non-determinism, an unintended code change, or a stale file."""
    path = Path(path)
    if path.exists() and path.read_text() != text:
        print(f"\n  !! REPRODUCIBILITY WARNING: {path.name} differs from the version already on disk "
              f"-- this run's output changed.\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_report(lines, filename):
    """Print a list-of-lines report and save it to results/<filename> (warns if it changed vs disk)."""
    text = "\n".join(lines)
    print(text)
    write_checked(RESULTS_DIR / filename, text + "\n")
    print(f"\nsaved -> {RESULTS_DIR / filename}")
