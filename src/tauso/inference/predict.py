"""Inference with the trained ASO-efficacy model.

Loads a trained XGBoost model by version and scores a DataFrame of features. The model predicts the
within-experiment efficacy *deviation* (the clean_exp target), so the output is a ranking score:
higher = more predicted knockdown relative to an experiment's mean. Rank candidates against each
other within a target; the absolute value is not a percent-inhibition prediction.

The boosters are large (~100 MB), so they are fetched from Zenodo on first use (like the feature
store) and cached under the data dir, rather than committed to git. Only the per-version feature
list ships in the package.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
import xgboost as xgb

from ..cli_utils import download_with_progress, verify_hash_or_exit
from ..data.data import get_data_dir

MODEL_DIR = Path(__file__).resolve().parent / "model"   # committed per-version feature lists
DEFAULT_VERSION = "v1"

# Zenodo record holding the trained boosters. Fill in the record id once the model is uploaded.
ZENODO_MODEL_RECORD = "21382921"
# Per version: the Zenodo model file + its md5 (pins the exact booster).
MODEL_FILES = {
    "v1": {"filename": "tauso_score_v1.json", "md5": "2a0b5bea756a7f623182f32d04eec3ad"},
}


def score_column(version=DEFAULT_VERSION):
    """Name of the score column this model version writes, e.g. 'tauso_score_v1'."""
    return f"tauso_score_{version}"


def ensure_model(version=DEFAULT_VERSION, force=False):
    """Path to the booster for `version`, fetching it from Zenodo into <data_dir>/models/ if absent
    (or if `force`). Verifies the md5 on every call so a truncated or wrong file fails loudly.

    `tauso setup-model` calls this to provision the booster ahead of time; `load_model` calls it
    on first use so inference works even without the setup step."""
    if version not in MODEL_FILES:
        raise KeyError(f"No model registered for version {version!r}.")
    spec = MODEL_FILES[version]
    dest = Path(get_data_dir()) / "models" / spec["filename"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    if force or not dest.exists():
        url = f"https://zenodo.org/api/records/{ZENODO_MODEL_RECORD}/files/{spec['filename']}/content"
        download_with_progress(url, str(dest), label=f"Downloading {spec['filename']}")
    verify_hash_or_exit(str(dest), spec["md5"], algo="md5")
    return dest


@lru_cache(maxsize=None)
def load_model(version=DEFAULT_VERSION):
    """Return (booster, feature_names) for `tauso_score_<version>`. The booster is fetched from
    Zenodo on first use (cached under the data dir); the feature list ships in the package."""
    booster = xgb.Booster()
    booster.load_model(str(ensure_model(version)))
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
