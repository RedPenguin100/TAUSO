import logging

import numpy as np
import pandas as pd
from notebooks.consts import OLIGO_CSV_PROCESSED_AVERAGED
from notebooks.features.feature_extraction import load_all_features

from tauso.data.consts import CANONICAL_GENE, CELL_LINE, INHIBITION, VOLUME
from tauso.timer import Timer

logger = logging.getLogger(__name__)


def get_data_path(version):
    if version == "oligo":
        return OLIGO_CSV_PROCESSED_AVERAGED
    raise ValueError(f"Invalid version: {version}")


def create_tauso_split(data: pd.DataFrame, val_frac: float = 0.15, test_frac: float = 0.15,
                       seed: int = 42) -> pd.DataFrame:
    """
    Create a stratified temporal split as an alternative to the OligoAI "split" column.

    Strategy: within each gene × cell-line cohort, sort by index_oligo (proxy for temporal
    order) and assign the first (1 - val_frac - test_frac) fraction to train, next val_frac
    to val, and final test_frac to test. Cohorts with < 3 rows go entirely to train.

    Returns the dataframe with a new "split" column (overwrites existing).
    """
    rng = np.random.default_rng(seed)
    train_frac = 1.0 - val_frac - test_frac
    split = pd.Series("train", index=data.index)

    for _, group in data.groupby([CANONICAL_GENE, CELL_LINE]):
        idx = group.sort_values("index_oligo").index
        n = len(idx)
        if n < 3:
            continue  # all train
        n_val  = max(1, round(n * val_frac))
        n_test = max(1, round(n * test_frac))
        n_train = n - n_val - n_test
        if n_train < 1:
            continue
        split.loc[idx[:n_train]]              = "train"
        split.loc[idx[n_train:n_train+n_val]] = "val"
        split.loc[idx[n_train+n_val:]]        = "test"

    data = data.copy()
    data["split"] = split

    counts = data["split"].value_counts()
    logger.info("Tauso split | train=%d val=%d test=%d",
                counts.get("train", 0), counts.get("val", 0), counts.get("test", 0))
    _ = rng  # seed kept for reproducibility but stratified split is deterministic
    return data


def _drop_legacy_other_column(final_data: pd.DataFrame) -> pd.DataFrame:
    """Backward-compat shim for caches built before populate_transfection started
    emitting only three columns with NaN-on-unknown semantics. Drops the legacy
    'Other' column if present and propagates NaN into Gymnosis / Lipofection /
    Electroporation for rows where Other == 1. Once the wide cache is regenerated
    this function becomes a no-op.
    """
    if "Other" not in final_data.columns:
        return final_data

    final_data = final_data.copy()
    other_mask = final_data["Other"].astype(bool)
    n_other = int(other_mask.sum())

    for col in ("Gymnosis", "Lipofection", "Electroporation"):
        if col in final_data.columns:
            final_data[col] = final_data[col].astype("float64")
            final_data.loc[other_mask, col] = np.nan

    final_data = final_data.drop(columns=["Other"])
    logger.info("Legacy 'Other' column dropped | %d rows had Other=1, now NaN on the three indicators", n_other)
    return final_data


def load_and_validate_final_data(version="oligo", load_competition=False, split_source="oligoai"):
    """
    Loads features and metadata, ensures shared columns are identical,
    and returns the merged DataFrame along with the final feature list.

    Transfection one-hots are always Electroporation / Gymnosis / Lipofection with
    NaN on rows whose transfection_method isn't one of those three (populate_transfection
    enforces this). Legacy caches that still carry an 'Other' column are normalized
    via _drop_legacy_other_column after the merge.
    """
    index = f"index_{version}"

    data_path = get_data_path(version)
    # 1. Load the data sources
    with Timer("Loading features"):
        loaded_features = load_all_features(version=version, load_competition=load_competition)
    data = pd.read_csv(data_path)

    # 2. Identify all columns that exist in both frames (other than the join key).
    common_cols = sorted(set(loaded_features.columns) & set(data.columns) - {index})

    # 3. Validation: ensure every shared column has identical values across sources.
    # Compare values not dtypes — the wide cache compresses some int columns to int8
    # for storage, while the same value in the source CSV stays int64.
    for col in common_cols:
        check_df = pd.merge(loaded_features[[index, col]], data[[index, col]], on=index)
        a, b = check_df[f"{col}_x"], check_df[f"{col}_y"]
        if a.isna().equals(b.isna()) and (a.dropna().to_numpy() == b.dropna().to_numpy()).all():
            del check_df
            continue
        raise ValueError(f"Integrity Check Failed: '{col}' differs between sources.")

    # 4. Drop the shared columns from `data` so the merge produces no _x/_y suffixes
    logger.info(f"Merge dedup | {len(common_cols)} shared columns dropped from data: {common_cols}")
    data_to_merge = data.drop(columns=common_cols)

    final_data = pd.merge(loaded_features, data_to_merge, on=index)

    # 5. Apply split
    if split_source == "tauso":
        logger.info("Creating tauso split (stratified temporal per gene×cell-line cohort)...")
        final_data = create_tauso_split(final_data)
    elif split_source != "oligoai":
        raise ValueError(f"Unknown split_source '{split_source}'. Use 'oligoai' or 'tauso'.")

    # 5b. Drop the legacy 'Other' column if a pre-NaN cache version sneaks in.
    final_data = _drop_legacy_other_column(final_data)

    # 6. Define Final Feature List
    features_to_ignore = [index, INHIBITION, "inhibition_percent", "dosage", "split"]

    features = [
        col for col in final_data.select_dtypes(include=["number"]).columns if col not in features_to_ignore
    ] + [VOLUME]

    # Ensure uniqueness and sort
    features = sorted(list(set(features)))

    return final_data, features


if __name__ == "__main__":
    from tauso.data.consts import *

    final_data, features = load_and_validate_final_data(version="oligo")
