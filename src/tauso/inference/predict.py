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
from typing import FrozenSet, List, Sequence, Tuple

import numpy as np
import pandas as pd
import xgboost as xgb

from ..data.data import get_data_dir

logger = logging.getLogger(__name__)

MODEL_DIR = Path(__file__).resolve().parent / "model"  # committed per-version feature lists
DEFAULT_VERSION = "v1"

# Zenodo record holding the trained boosters.
ZENODO_MODEL_RECORD = "22543712"
# Per version: the Zenodo model file + its md5 (pins the exact booster). `filename` must be the exact
# name of the file uploaded to the Zenodo record; it is fetched and cached under that same name.
MODEL_FILES = {
    "v1": {"filename": "tauso_score_v1.json", "md5": "7db62ebd7cfa08affd49b16a5b7c1938"},
}


def score_column(version: str = DEFAULT_VERSION) -> str:
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
def load_finite_features(version: str = DEFAULT_VERSION) -> FrozenSet[str]:
    """
    Load all features that can't be NaN, so we check validity.
    """
    path = MODEL_DIR / f"tauso_score_{version}.finite.txt"
    return frozenset(path.read_text().split())


Offenders = List[Tuple[str, int]]


def find_nonfinite(X: np.ndarray, features: Sequence[str], guarded: FrozenSet[str]) -> Offenders:
    """(feature, n_candidates) worst first. NaN counts only for `guarded`; infinity for all."""
    nonfinite = ~np.isfinite(X)
    if not nonfinite.any():
        return []
    is_guarded = np.array([name in guarded for name in features])
    counts = np.where(is_guarded, nonfinite, np.isinf(X)).sum(axis=0)
    found = [(name, int(n)) for name, n in zip(features, counts) if n]
    found.sort(key=lambda item: (-item[1], item[0]))
    return found


def classify_nonfinite(offenders: Offenders, n_candidates: int) -> Tuple[Offenders, Offenders]:
    """Split into (failures, gaps): failed for some candidates, vs missing for all of them."""
    failures = [item for item in offenders if item[1] < n_candidates]
    gaps = [item for item in offenders if item[1] == n_candidates]
    return failures, gaps


def _list_features(offenders: Offenders, limit: int = 5) -> str:
    shown = ", ".join(f"{name} ({count} candidates)" for name, count in offenders[:limit])
    return shown + (f" and {len(offenders) - limit} more" if len(offenders) > limit else "")


def check_scorable(X: np.ndarray, features: Sequence[str], version: str = DEFAULT_VERSION, strict: bool = True) -> None:
    """Raise if a guarded feature failed for part of the batch; warn if one is missing for all."""
    failures, gaps = classify_nonfinite(find_nonfinite(X, features, load_finite_features(version)), len(X))
    if gaps:
        logger.warning(
            f"{len(gaps)} of the model's {len(features)} features are non-finite for all {len(X)} "
            f"candidates: {_list_features(gaps)}. They are unavailable for this target rather than "
            f"broken, so treat the resulting order as approximate."
        )
    if failures:
        message = (
            f"{len(failures)} of the model's {len(features)} features failed to compute for part of "
            f"the {len(X)} candidates: {_list_features(failures)}. Those candidates would be ranked "
            f"on an arbitrary default. Re-run the feature pipeline for them or drop them; pass "
            f"strict=False to score them regardless."
        )
        if strict:
            raise ValueError(message)
        logger.warning(message)


def predict(features_df: pd.DataFrame, version: str = DEFAULT_VERSION, strict: bool = True) -> np.ndarray:
    """Score `features_df` -> ranking scores, higher = better predicted knockdown.

    Selects and orders to the model's feature list, so extra columns are ignored and a missing one
    is an error. See `check_scorable` for `strict`."""
    booster, features = load_model(version)
    missing = [f for f in features if f not in features_df.columns]
    if missing:
        raise ValueError(
            f"features_df is missing {len(missing)} of the model's {len(features)} features, e.g. {missing[:5]}"
        )
    X = features_df[features].to_numpy(np.float64)
    check_scorable(X, features, version, strict)
    return booster.predict(xgb.DMatrix(X, feature_names=features))


def score(features_df: pd.DataFrame, version: str = DEFAULT_VERSION, strict: bool = True) -> pd.DataFrame:
    """Return a copy of `features_df` with the model's score added as `tauso_score_<version>`, sorted best-first."""
    out = features_df.copy()
    col = score_column(version)
    out[col] = predict(out, version, strict=strict)
    return out.sort_values(col, ascending=False, kind="stable")
