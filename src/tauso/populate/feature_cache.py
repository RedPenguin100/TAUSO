import logging
import os

import numpy as np
import pandas as pd

from ..cli_utils import download_with_progress, echo_ok, verify_hash_or_exit
from ..data.data import get_data_dir

logger = logging.getLogger(__name__)

# Per run: the wide cache parquet on Zenodo + the index column it uses. Add a new entry
# (e.g. "oligo_v2") to ship a new version; runs can share `index_col` across versions.
ZENODO_FEATURES_RECORD = "20420643"
FEATURE_CACHE_FILES = {
    "oligo": {
        "filename": "oligo_features.parquet",
        "md5": "5978c7eae34aff7125d89839c16614e5",
        "index_col": "index_oligo",
    },
}


def _spec(run):
    if run not in FEATURE_CACHE_FILES:
        raise KeyError(f"No feature cache registered for run '{run}'.")
    return FEATURE_CACHE_FILES[run]


def feature_store_dir(run):
    d = os.path.join(get_data_dir(), "features", run)
    os.makedirs(d, exist_ok=True)
    return d


def cache_path(run):
    return os.path.join(feature_store_dir(run), _spec(run)["filename"])


def index_col(run):
    return _spec(run)["index_col"]


def ensure_cache(run, force=False):
    """Fetch the wide feature cache for `run` from Zenodo if not present locally."""
    dest = cache_path(run)
    expected = _spec(run)["md5"]
    if os.path.exists(dest) and not force:
        verify_hash_or_exit(dest, expected, algo="md5")
        return dest
    url = f"https://zenodo.org/api/records/{ZENODO_FEATURES_RECORD}/files/{_spec(run)['filename']}/content"
    download_with_progress(url, dest, label=f"Downloading {_spec(run)['filename']}")
    verify_hash_or_exit(dest, expected, algo="md5")
    echo_ok(f"Downloaded feature cache: {dest}")
    return dest


def save_feature_internal(df, feature_name, overwrite=False, version=None, saved_dir_func=None):
    if version is None:
        return  # version=None means "don't cache"
    if saved_dir_func is None:
        raise ValueError(f"saved_dir_func is required to cache feature '{feature_name}'.")

    feature_dir = saved_dir_func(version)
    os.makedirs(feature_dir, exist_ok=True)
    idx = _spec(version)["index_col"] if version in FEATURE_CACHE_FILES else f"index_{version}"
    sub_df = df[[idx, feature_name]].copy()
    file_path = os.path.join(feature_dir, f"{feature_name}.parquet")

    if os.path.exists(file_path):
        diff_rows, col_new, col_old = check_saved_features_conflicts(sub_df, file_path, feature_name, idx)
        if diff_rows.empty:
            existing = pd.read_parquet(file_path)
            new_rows = sub_df[~sub_df[idx].isin(existing[idx])]
            if not new_rows.empty:
                pd.concat([existing, new_rows], ignore_index=True).to_parquet(file_path, index=False)
                logger.info("Appended %d new rows to %s.parquet", len(new_rows), feature_name)
            return
        if overwrite:
            sub_df.to_parquet(file_path, index=False)
            logger.info("Overwrote feature: %s", feature_name)
            return
        raise ValueError(f"Collision detected for feature '{feature_name}'. Aborting.")

    sub_df.to_parquet(file_path, index=False)
    logger.info("Saved feature: %s", feature_name)


def check_saved_features_conflicts(new_sub_df, file_path, feature_name, index):
    existing_df = pd.read_parquet(file_path)
    comparison = pd.merge(new_sub_df, existing_df, on=index, suffixes=("_new", "_existing"), how="inner")
    col_new = f"{feature_name}_new"
    col_old = f"{feature_name}_existing"
    if pd.api.types.is_numeric_dtype(comparison[col_new]):
        diff_mask = ~np.isclose(comparison[col_new], comparison[col_old], rtol=1e-05, atol=1e-08, equal_nan=True)
    else:
        diff_mask = comparison[col_new].fillna("NA").astype(str) != comparison[col_old].fillna("NA").astype(str)
    return comparison[diff_mask], col_new, col_old
