import pytest
import pandas as pd
from notebooks.utils.oligo_ai_process import process_chemistry
from notebooks.competitors.oligo_ai.parse_chemistry import parse_sugar_string
from tauso.populate.populate_rnase import populate_rnase_features
from tauso.data.consts import CHEMICAL_PATTERN, MODIFICATION


@pytest.fixture(scope="module")
def chemistry_data(base_data):
    """Applies chemistry parsing to the universal base_data."""
    data = base_data.copy()
    data[CHEMICAL_PATTERN] = data['sugar_mods'].apply(parse_sugar_string)

    results = data['sugar_mods'].apply(process_chemistry)
    data[[CHEMICAL_PATTERN, MODIFICATION]] = pd.DataFrame(results.tolist(), index=data.index)

    return data


@pytest.fixture
def sampled_chem_data(request, chemistry_data):
    n_samples = getattr(request, 'param', 1000)
    actual_samples = min(n_samples, len(chemistry_data))
    return chemistry_data.sample(n=actual_samples, random_state=42).copy()


@pytest.mark.parametrize("sampled_chem_data", [10000], indirect=True)
def test_rnase_features_regression(sampled_chem_data, dataframe_regression):
    data = sampled_chem_data.copy()

    processed_data, rnase_features = populate_rnase_features(data)

    # Assuming rnase_features is a list of the new column names
    dataframe_regression.check(processed_data[rnase_features])