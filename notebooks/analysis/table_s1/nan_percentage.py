"""Table S1: how often each of the model's features is missing.

NaN is meaningful here -- nothing is imputed, and XGBoost learns a default branch -- so the
supplementary table has to state the missing fraction per feature. Writes out/nan_percentage.csv
with one row per model feature so the table can be diffed after a feature-pipeline change.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from notebooks.analysis.table_s1.consts import NAN_PERCENTAGE_CSV
from notebooks.models.common import load_dataset

from tauso.inference import DEFAULT_VERSION, get_model_feature_names

# The count Table S1 reports. It is a claim in the paper, so pin it rather than infer it.
EXPECTED_FEATURES = 485


def model_features(version: str = DEFAULT_VERSION) -> list[str]:
    """The shipped model's feature list, which is what Table S1 enumerates."""
    features = get_model_feature_names(version)
    if len(features) != EXPECTED_FEATURES or len(set(features)) != EXPECTED_FEATURES:
        raise AssertionError(
            f"model {version} lists {len(features)} features ({len(set(features))} unique), "
            f"expected exactly {EXPECTED_FEATURES}"
        )
    return features


def missingness(df: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    """Per-feature NaN count and percentage, sorted by name so the file diffs cleanly."""
    absent = [f for f in features if f not in df.columns]
    if absent:
        raise AssertionError(f"{len(absent)} model features missing from the dataset: {absent}")

    counts = df[features].isna().sum()
    table = pd.DataFrame(
        {
            "feature": counts.index,
            "n_missing": counts.to_numpy(),
            "pct_missing": (100 * counts / len(df)).round(2).to_numpy(),
        }
    ).sort_values("feature", ignore_index=True)

    if len(table) != EXPECTED_FEATURES or table["feature"].nunique() != EXPECTED_FEATURES:
        raise AssertionError(
            f"table has {len(table)} rows ({table['feature'].nunique()} distinct features), "
            f"expected exactly {EXPECTED_FEATURES}"
        )
    return table


def main():
    df, _ = load_dataset()
    table = missingness(df, model_features())

    NAN_PERCENTAGE_CSV.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(NAN_PERCENTAGE_CSV, index=False)

    ever_missing = table[table["n_missing"] > 0].sort_values(["n_missing", "feature"], ascending=[False, True])
    print(f"{len(df):,} rows, {len(table)} model features")
    print(f"{len(ever_missing)} ever missing, {len(table) - len(ever_missing)} never missing\n")
    print(ever_missing.to_string(index=False))
    print(f"\nall {len(table)} features -> {NAN_PERCENTAGE_CSV}")


if __name__ == "__main__":
    main()
