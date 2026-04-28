import os

import numpy as np
import pandas as pd


def save_feature_internal(df, feature_name, overwrite=False, version=None, saved_dir_func=None):
    if saved_dir_func is None:
        raise ValueError(
            f"Please provide a function for saving feature '{feature_name}' that accepts a version and outputs a directory."
        )

    if version is None:
        return

    feature_dir = saved_dir_func(version)

    os.makedirs(feature_dir, exist_ok=True)

    if version is None:
        index = "index"
    else:
        index = "index_" + version
    sub_df = df[[index, feature_name]].copy()
    file_path = os.path.join(feature_dir, f"{feature_name}.csv")

    file_exists = os.path.exists(file_path)

    if file_exists:
        diff_rows, col_new, col_old = check_saved_features_conflicts(sub_df, file_path, feature_name, index)

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


def check_saved_features_conflicts(new_sub_df, file_path, feature_name, index):
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
