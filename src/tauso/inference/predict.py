"""Inference with the trained ASO-efficacy model.

Loads a bundled XGBoost model by version and scores a DataFrame of features. The model predicts the
within-experiment efficacy *deviation* (the clean_exp target), so the output is a ranking score:
higher = more predicted knockdown relative to an experiment's mean. Rank candidates against each
other within a target; the absolute value is not a percent-inhibition prediction.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
import xgboost as xgb

MODEL_DIR = Path(__file__).resolve().parent / "model"
DEFAULT_VERSION = "v1"


def score_column(version=DEFAULT_VERSION):
    """Name of the score column this model version writes, e.g. 'tauso_score_v1'."""
    return f"tauso_score_{version}"


@lru_cache(maxsize=None)
def load_model(version=DEFAULT_VERSION):
    """Return (booster, feature_names) for the bundled `tauso_score_<version>` model. Cached per version."""
    stem = f"tauso_score_{version}"
    booster = xgb.Booster()
    booster.load_model(str(MODEL_DIR / f"{stem}.json"))
    features = (MODEL_DIR / f"{stem}.features.txt").read_text().splitlines()
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
