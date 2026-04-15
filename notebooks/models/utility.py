import pandas as pd

from notebooks.consts import OLIGO_CSV_PROCESSED_AVERAGED
from notebooks.features.feature_extraction import load_all_features
from tauso.data.consts import INHIBITION, VOLUME

def get_data_path(version):
    if version == 'oligo':
        return OLIGO_CSV_PROCESSED_AVERAGED
    raise ValueError(f"Invalid version: {version}")


def load_and_validate_final_data(version='oligo'):
    """
    Loads features and metadata, ensures shared columns are identical,
    and returns the merged DataFrame along with the final feature list.
    """
    index = f"index_{version}"

    data_path = get_data_path(version)
    # 1. Load the data sources
    loaded_features = load_all_features(version=version, n_jobs=-1)
    data = pd.read_csv(data_path)

    # 2. Identify common columns to validate (excluding the join key)
    # Based on your data: sense_start, sense_start_frokm_end, sense_length, etc.
    common_cols = [
        'sense_start', 'sense_start_from_end', 'sense_length',
        'sense_exon', 'sense_intron', 'sense_utr', 'sense_type'
    ]

    # 3. Validation: Ensure both have the columns and they are identical
    for col in common_cols:
        if col not in loaded_features.columns or col not in data.columns:
            raise KeyError(f"Column '{col}' missing from one of the dataframes.")

        # Merge temporarily on index to align rows for comparison
        check_df = pd.merge(
            loaded_features[[index, col]],
            data[[index, col]],
            on=index
        )

        if not check_df[f"{col}_x"].equals(check_df[f"{col}_y"]):
            raise ValueError(f"Integrity Check Failed: '{col}' differs between sources.")

    # 4. Merge: Drop redundant columns from 'data' first to avoid _x/_y
    data_to_merge = data.drop(columns=common_cols)
    final_data = pd.merge(loaded_features, data_to_merge, on=index)

    # 5. Define Final Feature List
    features_to_ignore = [index, INHIBITION, 'inhibition_percent', 'dosage', 'split']

    features = [
        col for col in final_data.select_dtypes(include=['number']).columns
        if col not in features_to_ignore
    ] + [VOLUME]

    # Ensure uniqueness and sort
    features = sorted(list(set(features)))

    return final_data, features