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


def _apply_other_as_nan(final_data: pd.DataFrame) -> pd.DataFrame:
    """
    Re-encode the 'Other' delivery method as missing-data semantics rather than
    a category. Where Other == 1:
      - Set Gymnosis / Lipofection / Electroporation to NaN
      - Add a single delivery_unknown flag (==1)
      - Drop the Other column
    For rows where Other == 0, delivery_unknown = 0 and the three indicators stay as-is.

    Rationale: 'Other' bundles >5 distinct real delivery methods (free uptake, LNP,
    in-vivo routes, etc). Treating it as one category lets the model overfit on
    spurious patterns. Treating it as missing data forces the model to rely on
    actual chemistry/expression features while keeping a single OOD flag.
    """
    if "Other" not in final_data.columns:
        logger.info("'Other' column not present; skipping NaN re-encoding.")
        return final_data

    final_data = final_data.copy()
    other_mask = final_data["Other"].astype(bool)
    n_other = int(other_mask.sum())

    for col in ("Gymnosis", "Lipofection", "Electroporation"):
        if col in final_data.columns:
            final_data[col] = final_data[col].astype("float64")
            final_data.loc[other_mask, col] = np.nan

    final_data["delivery_unknown"] = other_mask.astype(int)
    final_data = final_data.drop(columns=["Other"])

    logger.info("'Other' re-encoded as NaN | %d rows flagged as delivery_unknown=1", n_other)
    return final_data


def load_and_validate_final_data(version="oligo", load_competition=False, split_source="oligoai",
                                 other_encoding="onehot"):
    """
    Loads features and metadata, ensures shared columns are identical,
    and returns the merged DataFrame along with the final feature list.

    other_encoding:
      - 'onehot' (default): keep Other as a 1/0 column alongside the others
      - 'nan': when Other == 1, NaN out Gymnosis/Lipofection/Electroporation,
               add delivery_unknown flag, drop Other column
    """
    index = f"index_{version}"

    data_path = get_data_path(version)
    # 1. Load the data sources
    with Timer("Loading features"):
        loaded_features = load_all_features(version=version, load_competition=load_competition)
    data = pd.read_csv(data_path)

    # 2. Identify common columns to validate (excluding the join key)
    # Based on your data: sense_start, sense_start_from_end, sense_length, etc.
    common_cols = [
        "sense_start",
        "sense_start_from_end",
        "sense_length",
        "sense_exon",
        "sense_intron",
        "sense_utr",
        "sense_type",
    ]

    # 3. Validation: Ensure both have the columns and they are identical
    for col in common_cols:
        if col not in loaded_features.columns or col not in data.columns:
            raise KeyError(f"Column '{col}' missing from one of the dataframes.")

        # Merge temporarily on index to align rows for comparison
        check_df = pd.merge(loaded_features[[index, col]], data[[index, col]], on=index)

        if not check_df[f"{col}_x"].equals(check_df[f"{col}_y"]):
            raise ValueError(f"Integrity Check Failed: '{col}' differs between sources.")
        del check_df

    # 4. Merge: Drop redundant columns from 'data' first to avoid _x/_y
    data_to_merge = data.drop(columns=common_cols)

    final_data = pd.merge(loaded_features, data_to_merge, on=index)

    # 5. Apply split
    if split_source == "tauso":
        logger.info("Creating tauso split (stratified temporal per gene×cell-line cohort)...")
        final_data = create_tauso_split(final_data)
    elif split_source != "oligoai":
        raise ValueError(f"Unknown split_source '{split_source}'. Use 'oligoai' or 'tauso'.")

    # 5b. Re-encode 'Other' delivery if requested
    if other_encoding == "nan":
        final_data = _apply_other_as_nan(final_data)
    elif other_encoding != "onehot":
        raise ValueError(f"Unknown other_encoding '{other_encoding}'. Use 'onehot' or 'nan'.")

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
