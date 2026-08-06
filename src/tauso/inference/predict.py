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


def predict(features_df, version=DEFAULT_VERSION):
    """Score `features_df` with model `version` -> np.ndarray of efficacy ranking scores (higher =
    better predicted knockdown). The frame must contain the model's feature columns; they are
    selected and ordered to the model's feature list, so extra columns are ignored and a missing
    feature is an error."""
    booster, features = load_model(version)
    missing = [f for f in features if f not in features_df.columns]
    if missing:
        raise ValueError(
            f"features_df is missing {len(missing)} of the model's {len(features)} features, e.g. {missing[:5]}"
        )
    X = features_df[features].to_numpy(np.float64)
    return booster.predict(xgb.DMatrix(X, feature_names=features))


def score(features_df, version=DEFAULT_VERSION):
    """Return a copy of `features_df` with the model's score added as `tauso_score_<version>`, sorted best-first."""
    out = features_df.copy()
    col = score_column(version)
    out[col] = predict(out, version)
    return out.sort_values(col, ascending=False, kind="stable")
