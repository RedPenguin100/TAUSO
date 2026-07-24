"""Inference with the trained ASO-efficacy model.

Loads a trained XGBoost model by version and scores a DataFrame of features. The model predicts the
within-experiment efficacy *deviation* (the clean_exp target), so the output is a ranking score:
higher = more predicted knockdown relative to an experiment's mean. Rank candidates against each
other within a target; the absolute value is not a percent-inhibition prediction.

The model is large (~100 MB) so it is not committed to git; download it with `tauso setup-model`.
"""

import logging
from functools import lru_cache
from pathlib import Path

import numpy as np
import xgboost as xgb

from ..data.data import get_data_dir

logger = logging.getLogger(__name__)

MODEL_DIR = Path(__file__).resolve().parent / "model"  # committed per-version feature lists
DEFAULT_VERSION = "v1"

# Zenodo record holding the trained boosters. Fill in the record id once the model is uploaded.
ZENODO_MODEL_RECORD = "21394367"
# Per version: the Zenodo model file + its md5 (pins the exact booster). `filename` must be the exact
# name of the file uploaded to the Zenodo record; it is fetched and cached under that same name.
MODEL_FILES = {
    "v1": {"filename": "deploy_model_best_clean_exp_trainval_seed1.json", "md5": "65e184e0111b5d4365d1b7d608ae9024"},
}


def score_column(version=DEFAULT_VERSION):
    """Name of the score column this model version writes, e.g. 'tauso_score_v1'."""
    return f"tauso_score_{version}"


def model_path(version=DEFAULT_VERSION):
    """Local path to the booster for `version`, under <data_dir>/models/. Does not download it;
    run `tauso setup-model` to provision it."""
    if version not in MODEL_FILES:
        raise KeyError(f"No model registered for version {version!r}.")
    return Path(get_data_dir()) / "models" / MODEL_FILES[version]["filename"]


@lru_cache(maxsize=None)
def load_model(version=DEFAULT_VERSION):
    """Return (booster, feature_names) for `tauso_score_<version>`. The booster must already be
    present locally (run `tauso setup-model` to download it); the feature list ships in the package."""
    path = model_path(version)
    if not path.exists():
        message = f"Model '{version}' not found at {path}. Run `tauso setup-model` to download it."
        logger.error(message)
        raise FileNotFoundError(message)
    booster = xgb.Booster()
    booster.load_model(str(path))
    features = (MODEL_DIR / f"tauso_score_{version}.features.txt").read_text().splitlines()
    return booster, features


@lru_cache(maxsize=None)
def load_finite_features(version=DEFAULT_VERSION):
    """The model features that are never missing in the training data, as a frozenset.

    Part of the feature set is conditionally defined -- codon scores exist only for windows inside a
    CDS, cEt hybridization terms only for cEt chemistries, splice-junction distances only near a
    junction -- and the booster learned an explicit missing-value branch for each. A NaN there is a
    real input the model handles.

    For the features named here it never saw a missing value, so a NaN is either a failed feature
    computation (a fold, RIsearch or hybridization call that returned nothing for one candidate) or an
    annotation that does not resolve for the target at all. XGBoost would still route it down an
    arbitrary default branch and return a plausible score, which is why `predict` flags it."""
    path = MODEL_DIR / f"tauso_score_{version}.finite.txt"
    return frozenset(f for f in path.read_text().splitlines() if f.strip())


def _find_nonfinite(X, features, version):
    """Non-finite values in `X` the model cannot interpret, as [(feature, n_candidates), ...] worst
    first. NaN counts only for the always-finite features; +-inf counts for all of them (no feature
    is ever infinite in training, so an infinity is out of distribution whatever the feature)."""
    nonfinite = ~np.isfinite(X)
    if not nonfinite.any():
        return []
    guarded = load_finite_features(version)
    infinite = np.isinf(X)
    found = []
    for col, name in enumerate(features):
        bad = nonfinite[:, col] if name in guarded else infinite[:, col]
        count = int(bad.sum())
        if count:
            found.append((name, count))
    found.sort(key=lambda item: -item[1])
    return found


def _classify_nonfinite(offenders, n_candidates):
    """Split `_find_nonfinite` output into (failures, gaps).

    A feature non-finite for only *some* candidates varies within one batch, which a feature that
    computed correctly cannot do: it failed for those candidates specifically, and they are the ones
    that get mis-scored relative to their neighbours. A feature non-finite for *every* candidate is
    instead unavailable for the target as a whole -- the gene-level annotations (ribosome profiling,
    mRNA half-life) have no value for a target outside the reference annotation, such as a custom
    `gene_sequence`. That perturbs the whole ranking rather than singling candidates out.

    With a single candidate the two are indistinguishable and everything reads as a gap."""
    failures = [item for item in offenders if item[1] < n_candidates]
    gaps = [item for item in offenders if item[1] == n_candidates]
    return failures, gaps


def _describe(offenders, limit=5):
    shown = ", ".join(f"{name} ({count} candidates)" for name, count in offenders[:limit])
    return shown + (f" and {len(offenders) - limit} more" if len(offenders) > limit else "")


def predict(features_df, version=DEFAULT_VERSION, strict=True):
    """Score `features_df` with model `version` -> np.ndarray of efficacy ranking scores (higher =
    better predicted knockdown). The frame must contain the model's feature columns; they are
    selected and ordered to the model's feature list, so extra columns are ignored and a missing
    feature is an error.

    A feature the model was never trained to see missing (see `load_finite_features`) that is
    non-finite for only part of the batch failed to compute for those candidates, and scoring them
    would return a plausible but meaningless number; that raises a ValueError, which `strict=False`
    downgrades to a warning. One that is non-finite for the whole batch is unavailable for the target
    rather than broken, and warns: the scores stay comparable to each other but the ranking is
    perturbed by an input the model has no branch for."""
    booster, features = load_model(version)
    missing = [f for f in features if f not in features_df.columns]
    if missing:
        raise ValueError(
            f"features_df is missing {len(missing)} of the model's {len(features)} features, e.g. {missing[:5]}"
        )
    X = features_df[features].to_numpy(np.float64)
    failures, gaps = _classify_nonfinite(_find_nonfinite(X, features, version), len(X))
    if gaps:
        logger.warning(
            f"{len(gaps)} of the model's {len(features)} features are non-finite for all {len(X)} "
            f"candidates: {_describe(gaps)}. They are unavailable for this target rather than broken -- "
            f"gene-level annotations resolve only for a target in the reference annotation, so a custom "
            f"gene_sequence has none. The model has no missing-value branch for them, so treat the "
            f"resulting order as approximate."
        )
    if failures:
        message = (
            f"{len(failures)} of the model's {len(features)} features failed to compute for part of the "
            f"{len(X)} candidates: {_describe(failures)}. The model has no missing-value branch for "
            f"these, so those candidates would be scored and ranked on an arbitrary default. Re-run the "
            f"feature pipeline for them or drop them; pass strict=False to score them regardless."
        )
        if strict:
            raise ValueError(message)
        logger.warning(message)
    return booster.predict(xgb.DMatrix(X, feature_names=features))


def score(features_df, version=DEFAULT_VERSION, strict=True):
    """Return a copy of `features_df` with the model's score added as `tauso_score_<version>`, sorted best-first."""
    out = features_df.copy()
    col = score_column(version)
    out[col] = predict(out, version, strict=strict)
    return out.sort_values(col, ascending=False, kind="stable")
