"""OligoGym's classical baselines on TAUSO's frozen split, scored by TAUSO's metrics.

Fits each (featurizer x model) pair from OligoGym on the train+val rows of the frozen OligoAI
split and scores the held-out test rows with `notebooks/models/evaluate.py`, so the numbers sit
in the same table as TAUSO's own. The models see only the oligo itself -- sequence and chemistry,
as HELM -- which is the point: it is the bar that TAUSO's genome-derived features have to clear.

Two targets are fitted, mirroring `notebooks/models/common.VARIANTS`:
  reg        raw inhibition_percent
  clean_exp  the deviation from each experiment's (custom_id) mean, which is what TAUSO ships

Run:
  python notebooks/competitors/oligogym/run_baselines.py
"""

import argparse
import inspect
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path[:0] = [str(HERE), str(REPO_ROOT), str(REPO_ROOT / "src")]

from helm_bridge import add_helm_column  # noqa: E402
from notebooks.models import common  # noqa: E402

from tauso.data.consts import (ASO_SEQUENCE, CANONICAL_GENE_NAME,  # noqa: E402
                               CELL_LINE, VOLUME_NM)

RESULTS_DIR = HERE / "results"


# Experiment covariates that can be appended to a featurizer's output. A HELM string describes
# the oligo and nothing about the experiment it was run in, so these are what the baselines are
# otherwise missing relative to the TAUSO model. Note they cannot move exp_med: all of them are
# constant within a custom_id, and exp_med ranks within one. They move gxc_* and globP.
EXTRA_GROUPS = {
    "dose": [VOLUME_NM, "density_cells_per_well"],
    "transfection": ["transfection_electroporation", "transfection_gymnosis",
                     "transfection_lipofection"],
    "cell_line": [CELL_LINE],
    "gene": [CANONICAL_GENE_NAME],
}
LOG_SCALED = {VOLUME_NM, "density_cells_per_well"}  # both span several orders of magnitude


def build_extras(trainval, test, groups):
    """Covariate columns for both sides, encoded on train+val only.

    Numeric columns are log-scaled where they span orders of magnitude and filled with the
    train+val median; categoricals become one-hots over the train+val categories, so a category
    seen only at test time encodes as all-zero rather than widening the matrix.
    """
    left, right, names = [], [], []
    for group in groups:
        for column in EXTRA_GROUPS[group]:
            if trainval[column].dtype.kind in "ifb":
                tr = trainval[column].to_numpy(np.float64)
                te = test[column].to_numpy(np.float64)
                if column in LOG_SCALED:
                    tr, te = np.log1p(tr), np.log1p(te)
                fill = np.nanmedian(tr)
                left.append(np.nan_to_num(tr, nan=fill)[:, None])
                right.append(np.nan_to_num(te, nan=fill)[:, None])
                names.append(column)
            else:
                categories = sorted(trainval[column].dropna().unique())
                left.append(pd.get_dummies(trainval[column]).reindex(
                    columns=categories, fill_value=0).to_numpy(np.float64))
                right.append(pd.get_dummies(test[column]).reindex(
                    columns=categories, fill_value=0).to_numpy(np.float64))
                names += [f"{column}={c}" for c in categories]
    return np.hstack(left), np.hstack(right), names


def build_featurizers(names):
    from oligogym.features import KMersCounts, OneHotEncoder, Thermodynamics

    available = {
        "kmers": lambda: KMersCounts(k=[1, 2, 3], modification_abundance=True),
        "onehot": lambda: OneHotEncoder(),
        "thermo": lambda: Thermodynamics(),
    }
    return {n: available[n]() for n in names}


def build_models(names, seed):
    from oligogym.models import (LinearModel, NearestNeighborsModel, RandomForestModel,
                                 XGBoostModel)

    available = {
        "linear": lambda: LinearModel(task="regression", type="ridge"),
        "knn": lambda: NearestNeighborsModel(task="regression"),
        "rf": lambda: RandomForestModel(task="regression", n_jobs=-1, random_state=seed),
        "xgboost": lambda: XGBoostModel(task="regression", n_jobs=-1, random_state=seed),
    }
    return {n: available[n] for n in names}


def as_array(features):
    """Featurizers return either an ndarray or a DataFrame; models want a 2-D float array."""
    if isinstance(features, pd.DataFrame):
        features = features.to_numpy()
    features = np.asarray(features, dtype=np.float64)
    if features.ndim > 2:
        features = features.reshape(features.shape[0], -1)
    return features


def featurize(featurizer, trainval, test, pad_length):
    """Fit the featurizer on train+val only, then transform both sides.

    `Thermodynamics` is stateless -- its `transform` re-runs `fit_transform`, so its width follows
    the longest oligo in whichever list it is handed. It is given an explicit `pad_length` so both
    sides pad to the same width, and DataFrame output is reindexed onto the train columns so a
    column that survives `dropna` on one side but not the other cannot shift the matrix.
    """
    kwargs = {}
    if "pad_length" in inspect.signature(featurizer.fit_transform).parameters:
        kwargs["pad_length"] = pad_length

    x_trainval = featurizer.fit_transform(trainval["helm"].tolist(), **kwargs)
    x_test = featurizer.transform(test["helm"].tolist(), **kwargs)
    if isinstance(x_trainval, pd.DataFrame):
        x_test = x_test.reindex(columns=x_trainval.columns, fill_value=0)
    return as_array(x_trainval), as_array(x_test)


def save(record, tag):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / f"oligogym_baselines{tag}.json").write_text(json.dumps(record, indent=2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--featurizers", nargs="+", default=["kmers", "onehot", "thermo"])
    ap.add_argument("--models", nargs="+", default=["linear", "knn", "rf", "xgboost"])
    ap.add_argument("--variants", nargs="+", default=["reg", "clean_exp"])
    ap.add_argument("--extra", nargs="*", default=[], choices=sorted(EXTRA_GROUPS),
                    help="experiment covariates to append to every featurizer's output")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tag", default="", help="suffix for the output filenames")
    args = ap.parse_args()

    # The repo's dataset, frozen split and cohort_id, plus a HELM column for the featurizers.
    # The feature columns load_dataset returns are ignored: these baselines read only the oligo.
    df, _ = common.load_dataset()
    df = add_helm_column(df)
    trainval, test = common.split(df)
    print(f"{len(df)} rows | train+val {len(trainval)} | test {len(test)} | "
          f"{test['custom_id'].nunique()} test experiments, {test['cohort_id'].nunique()} cohorts", flush=True)

    pad_length = int(df[ASO_SEQUENCE].str.len().max())
    featurizers = build_featurizers(args.featurizers)
    models = build_models(args.models, args.seed)

    rows, record = {}, []
    for feature_name, featurizer in featurizers.items():
        t0 = time.time()
        x_trainval, x_test = featurize(featurizer, trainval, test, pad_length)
        if args.extra:
            extra_trainval, extra_test, extra_names = build_extras(trainval, test, args.extra)
            x_trainval = np.hstack([x_trainval, extra_trainval])
            x_test = np.hstack([x_test, extra_test])
        print(f"\n[{feature_name}] {x_trainval.shape[1]} features in {time.time() - t0:.0f}s"
              + (f" (+{len(extra_names)} covariates)" if args.extra else ""), flush=True)

        for variant in args.variants:
            y_trainval = common._target(trainval, common.VARIANTS[variant])
            for model_name, make_model in models.items():
                label = f"{feature_name}/{model_name}/{variant}"
                t0 = time.time()
                try:
                    model = make_model()
                    model.fit(x_trainval, y_trainval)
                    predictions = np.asarray(model.predict(x_test), dtype=np.float64).ravel()
                except Exception as exc:  # one dead combination must not lose the whole run
                    print(f"  {label:34s} FAILED: {type(exc).__name__}: {exc}", flush=True)
                    record.append({"featurizer": feature_name, "model": model_name,
                                   "variant": variant, "error": f"{type(exc).__name__}: {exc}"})
                    continue
                scores = common.metrics_from(predictions, test)
                rows[label] = scores
                record.append({"featurizer": feature_name, "model": model_name, "variant": variant,
                               "extra": list(args.extra),
                               "n_features": int(x_trainval.shape[1]), "seconds": round(time.time() - t0, 1),
                               **{k: float(v) for k, v in scores.items()}})
                print(f"  {label:34s} exp_med {scores['exp_med']:.4f}  gxc_med {scores['gxc_med']:.4f}  "
                      f"({time.time() - t0:.0f}s)", flush=True)
                save(record, args.tag)  # banked as we go, so a killed run keeps what it finished

    header = (f"OligoGym classical baselines on TAUSO's frozen split -- held-out test"
              f"   train+val={len(trainval)}   test={len(test)}"
              + (f"   + covariates: {', '.join(args.extra)}" if args.extra else ""))
    lines = [header, "", *common.metric_table("HELD-OUT TEST (sequence + chemistry only)", rows)]
    text = "\n".join(lines)
    print("\n" + text)

    save(record, args.tag)
    stem = f"oligogym_baselines{args.tag}"
    (RESULTS_DIR / f"{stem}.txt").write_text(text + "\n")
    print(f"\nsaved -> {RESULTS_DIR / stem}.txt / .json")


if __name__ == "__main__":
    main()
