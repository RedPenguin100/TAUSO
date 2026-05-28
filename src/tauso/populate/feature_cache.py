import logging
import os

import numpy as np
import pandas as pd

from ..cli_utils import download_with_progress, echo_ok, verify_hash_or_exit
from ..data.data import get_data_dir

logger = logging.getLogger(__name__)

# --- Feature store layout -------------------------------------------------------------
# A run's store lives at TAUSO_DATA_DIR/features/<run>/ and holds:
#   <run>_features.parquet : the wide, frozen feature CACHE (downloaded from Zenodo)
#   <feature>.parquet      : loose per-feature DEV shards (local; override/extend the cache)
# The cache is just a big multi-column parquet; load_all_features reads it plus every shard
# and outer-joins them (shards win on name clash), so adding/removing a feature is just
# writing/deleting one shard -- the cache is never rewritten in place.
ZENODO_FEATURES_RECORD = "20420643"  # https://zenodo.org/records/20420643
FEATURE_CACHE_FILES = {
    # run -> wide cache parquet in the Zenodo record + its md5 (both printed by pack_features.py)
    "oligo": {"filename": "oligo_features.parquet", "md5": "5978c7eae34aff7125d89839c16614e5"},
}


def feature_store_dir(run):
    d = os.path.join(get_data_dir(), "features", run)
    os.makedirs(d, exist_ok=True)
    return d


def cache_filename(run):
    spec = FEATURE_CACHE_FILES.get(run)
    return spec["filename"] if spec else f"{run}_features.parquet"


def cache_path(run):
    return os.path.join(feature_store_dir(run), cache_filename(run))


def ensure_cache(run, force=False):
    """Ensure the wide feature cache for `run` is present locally, fetching from Zenodo if not.

    Returns the local cache path. Raises FileNotFoundError if the run is unpublished
    (placeholder Zenodo coords) and no local copy exists.
    """
    dest = cache_path(run)
    spec = FEATURE_CACHE_FILES.get(run)

    if os.path.exists(dest) and not force:
        if spec and not spec["md5"].startswith("PLACEHOLDER"):
            verify_hash_or_exit(dest, spec["md5"], algo="md5")
        return dest

    if spec is None:
        raise FileNotFoundError(
            f"No feature cache registered for run '{run}'. Known runs: {', '.join(FEATURE_CACHE_FILES)}."
        )
    if ZENODO_FEATURES_RECORD.startswith("PLACEHOLDER") or spec["md5"].startswith("PLACEHOLDER"):
        raise FileNotFoundError(
            f"Feature cache '{run}' is not published yet (placeholder Zenodo coords).\n"
            f"  Build it with `python -m notebooks.features.pack_features --run {run}`, upload to Zenodo,\n"
            f"  and set the record id + md5 here -- or drop the parquet at:\n    {dest}"
        )

    url = f"https://zenodo.org/api/records/{ZENODO_FEATURES_RECORD}/files/{spec['filename']}/content"
    download_with_progress(url, dest, label=f"Downloading {spec['filename']}")
    verify_hash_or_exit(dest, spec["md5"], algo="md5")
    echo_ok(f"Downloaded feature cache: {dest}")
    return dest


def save_feature_internal(df, feature_name, overwrite=False, version=None, saved_dir_func=None):
    # version=None means "don't cache" (e.g. computing features for arbitrary new inputs).
    if version is None:
        return
    if saved_dir_func is None:
        raise ValueError(f"saved_dir_func is required to cache feature '{feature_name}'.")

    feature_dir = saved_dir_func(version)
    os.makedirs(feature_dir, exist_ok=True)

    index = "index_" + version
    sub_df = df[[index, feature_name]].copy()
    file_path = os.path.join(feature_dir, f"{feature_name}.parquet")

    if os.path.exists(file_path):
        diff_rows, col_new, col_old = check_saved_features_conflicts(sub_df, file_path, feature_name, index)

        if diff_rows.empty:
            # Values agree on overlapping rows: append any genuinely new index rows.
            existing = pd.read_parquet(file_path)
            new_rows = sub_df[~sub_df[index].isin(existing[index])]
            if not new_rows.empty:
                pd.concat([existing, new_rows], ignore_index=True).to_parquet(file_path, index=False)
                logger.info("File exists for '%s'. Added %d new rows.", feature_name, len(new_rows))
            else:
                logger.info("File exists for '%s' and values match (within tolerance). No action taken.", feature_name)
            return

        logger.warning("Conflict detected for feature: %s", feature_name)
        logger.warning(
            "Found %d differing rows. Top 10 differences:\n%s",
            len(diff_rows),
            diff_rows[[index, col_new, col_old]].head(10).to_string(),
        )

        if overwrite:
            sub_df.to_parquet(file_path, index=False)
            logger.info("Overwrote feature: %s", feature_name)
            return
        raise ValueError(f"Collision detected for feature '{feature_name}'. Aborting.")

    sub_df.to_parquet(file_path, index=False)
    logger.info("Saved feature: %s", feature_name)


def check_saved_features_conflicts(new_sub_df, file_path, feature_name, index):
    existing_df = pd.read_parquet(file_path)

    # Merge on the index column to compare aligned rows.
    comparison = pd.merge(new_sub_df, existing_df, on=index, suffixes=("_new", "_existing"), how="inner")

    col_new = f"{feature_name}_new"
    col_old = f"{feature_name}_existing"

    if pd.api.types.is_numeric_dtype(comparison[col_new]):
        is_close = np.isclose(comparison[col_new], comparison[col_old], rtol=1e-05, atol=1e-08, equal_nan=True)
        diff_mask = ~is_close
    else:
        # Normalize NaN/strings so the string "NA" doesn't collide with a real NaN.
        new_normalized = comparison[col_new].fillna("NA").astype(str)
        old_normalized = comparison[col_old].fillna("NA").astype(str)
        diff_mask = new_normalized != old_normalized

    diff_rows = comparison[diff_mask]
    return diff_rows, col_new, col_old
