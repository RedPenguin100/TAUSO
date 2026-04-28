import pytest

from tauso.populate.populate_fold import (
    populate_mfe_features,
    populate_sense_accessibility_batch,
)


@pytest.fixture
def mini_sampled_data(request, final_data):
    """Samples the fully processed DataFrame right before the test runs."""
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(final_data))
    return final_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_access(mini_sampled_data, gene_to_data, dataframe_regression):
    configurations = [
        {"flank": 120, "access": 20, "seeds": [13]},
        {"flank": 120, "access": 13, "seeds": [4, 6, 8]},
        {"flank": 120, "access": 20, "seeds": [4, 6, 8]},
        {"flank": 120, "access": 13, "seeds": [13, 26, 39]},
    ]

    feature_names = []
    for config in configurations:
        c_flank = config["flank"]
        c_access = config["access"]
        c_seeds = config["seeds"]

        print(f"Running: Flank={c_flank}, Access={c_access}, Seeds={c_seeds})...")

        mini_sampled_data, feature_name = populate_sense_accessibility_batch(
            mini_sampled_data,
            gene_to_data,
            batch_size=10,
            flank_size=c_flank,
            access_size=c_access,
            seed_sizes=c_seeds,
            n_jobs=32,
        )
        feature_names.append(feature_name)

    dataframe_regression.check(mini_sampled_data[["index_oligo"] + feature_names])


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_mfe(mini_sampled_data, gene_to_data, dataframe_regression):
    mini_sampled_data, feature_names = populate_mfe_features(mini_sampled_data, gene_to_data)

    dataframe_regression.check(mini_sampled_data[["index_oligo"] + feature_names])
