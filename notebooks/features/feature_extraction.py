import os

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from notebooks.consts import SAVED_FEATURES


def _get_saved_features_dir(version):
    if version is None:
        return SAVED_FEATURES
    elif version in ["v2", "oligo"]:
        return SAVED_FEATURES.with_name(f"{SAVED_FEATURES.name}_{version}")
    else:
        raise ValueError(f"Unknown version '{version}'")


COMPETITION = [
    "PFRED_PLS",
    "PFRED_SVM",
    "OW_Overall",
    "OW_Tm",
    "OW_Intra_Oligo",
    "OW_Duplex",
    "sfold_accessibility",
    "miranda_score",
    "miranda_energy",
    "oligo_ai_score",
]


def load_all_features(filenames=None, light=True, verbose=False, version=None, n_jobs=1, load_competition=False):
    feature_dir = _get_saved_features_dir(version)
    if not filenames:
        filenames = [f for f in os.listdir(feature_dir) if f.endswith(".csv") and not f.startswith(".")]
        filenames.sort()

    if not filenames:
        raise FileNotFoundError(f"No CSV files found in {feature_dir}")

    if light:
        if "Smiles.csv" in filenames:
            filenames.remove("Smiles.csv")

    if not load_competition:
        filenames = [f for f in filenames if f.removesuffix(".csv") not in COMPETITION]

    if verbose:
        print(f"Loading features from: {filenames}")

    if version is None:
        index = "index"
    else:
        index = "index_" + version

    if n_jobs == 1:
        dfs = [pd.read_csv(os.path.join(feature_dir, f), index_col=index) for f in filenames]
    else:
        dfs = Parallel(n_jobs=n_jobs, backend="threading")(
            delayed(pd.read_csv)(os.path.join(feature_dir, f), index_col=index) for f in filenames
        )

    merged_df = pd.concat(dfs, axis=1, join="outer").reset_index()

    return merged_df


def _check_conflicts(new_sub_df, file_path, feature_name, index):
    existing_df = pd.read_csv(file_path)

    # Merge on 'index' to compare aligned rows
    comparison = pd.merge(new_sub_df, existing_df, on=index, suffixes=("_new", "_existing"), how="inner")

    col_new = f"{feature_name}_new"
    col_old = f"{feature_name}_existing"

    # Check if the column is numeric (float/int)
    if pd.api.types.is_numeric_dtype(comparison[col_new]):
        is_close = np.isclose(comparison[col_new], comparison[col_old], rtol=1e-05, atol=1e-08, equal_nan=True)
        diff_mask = ~is_close
    else:
        # For strings/other types, fill actual NaNs with "NA" and ensure both are strings
        # This prevents the string "NA" from conflicting with the float NaN
        new_normalized = comparison[col_new].fillna("NA").astype(str)
        old_normalized = comparison[col_old].fillna("NA").astype(str)

        diff_mask = new_normalized != old_normalized

    diff_rows = comparison[diff_mask]
    return diff_rows, col_new, col_old


def save_feature(df, feature_name, overwrite=False, version=None):
    if version is None:
        return

    feature_dir = _get_saved_features_dir(version)

    os.makedirs(feature_dir, exist_ok=True)

    if version is None:
        index = "index"
    else:
        index = "index_" + version
    sub_df = df[[index, feature_name]].copy()
    file_path = os.path.join(feature_dir, f"{feature_name}.csv")

    file_exists = os.path.exists(file_path)

    if file_exists:
        diff_rows, col_new, col_old = _check_conflicts(sub_df, file_path, feature_name, index)

        if diff_rows.empty:
            # 1. Check for new indexes and append them efficiently
            existing_idx = pd.read_csv(file_path, usecols=[index])[index]
            new_rows = sub_df[~sub_df[index].isin(existing_idx)]

            if not new_rows.empty:
                new_rows.to_csv(file_path, mode="a", header=False, index=False)
                print(f"File exists for '{feature_name}'. Appended {len(new_rows)} new rows.")
            else:
                print(f"File exists for '{feature_name}' but values are identical (within tolerance). No action taken.")
            return

        print(f"!!! Conflict detected for feature: {feature_name} !!!")
        print(f"Found {len(diff_rows)} differing rows. Top 10 differences:")
        print(diff_rows[[index, col_new, col_old]].head(10))
        print("-" * 30)

        if overwrite:
            sub_df.to_csv(file_path, index=False)
            print(f"Overwrote feature: {feature_name}")
            return  # Prevent fall-through to the non-existent file block below
        else:
            # 2. Abort on collision
            raise ValueError(f"Collision detected for feature '{feature_name}'. Aborting.")

    # 3. Only reached if the file does not exist initially
    sub_df.to_csv(file_path, index=False)
    print(f"Saved feature: {feature_name}")


# TODO: fix this function
# def read_base_df():
#     all_data = pd.read_csv(DATA_PATH_NEW / 'data_asoptimizer_updated.11.7.csv',
#                            dtype={'index': int, 'ISIS': int, 'Target_gene': str,
#                                   'Cell_line': str, 'Density(cells_per_well)': float,
#                                   'Transfection': str, 'ASO_volume(nm)': float, 'Treatment_Period(hours)': float,
#                                   'Primer_probe_set': str, 'Sequence': str, 'Modification': str, 'Location': str,
#                                   'Chemical_Pattern': str, 'Linkage': str, 'Linage_Location' : str, 'Smiles' : str,
#                                   'Inhibition(%)' : float, 'seq_length' : int, 'Canonical Gene Name' : str,
#                                   'Cell line organism' : str, 'Transcript' : str, 'Location_in_sequence' : int,
#                                   'Location_div_by_length' : float, 'true_length_of_seq' : int, 'mod_scan' : int,
#                                   'cell_line_uniform' : str
#                                   })
#     return all_data
