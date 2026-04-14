import os
from typing import Dict

import pandas as pd
import numpy as np

import xgboost as xgb

import pytest

from notebooks.consts import NOTEBOOK_PATH
from tauso.aso_generation import generate_aso_features, generate_stub_data, Transfection
from tauso.genome.LocusInfo import LocusInfo
from tauso.populate.calculators.cache import AssetCache


@pytest.fixture(scope="session")
def cache():
    cache = AssetCache()
    # with Timer():
    #     cache.preload()
    return cache


TMSB10_PREMRNA = 'GTTTCTTGCTGCAGCAACGCGAGTGGGAGCACCAGGATCTCGGGCTCGGAACGAGACTGCACGGTGAGTGCGGCGCCGGGGCGGGGGGCCCACCCAGGGTGTGGTCGGATCCGGTGCACCGGGCGGGCGCGCCGCAACCGCGACAGGCGCCCTTCTCGGACCGGACGCAGGGGCCGGCGACCACGCCCTGGGACCGAGAAGAGGGGTGCGGGACGCGCCCAGATCCTCGGCCTTGGGGCTGCTCGGCACGCCTTGGCGCGAGTGCCACGTCGAGAGGCGTCGGCGGGGAGCGCGGAAGGGGACGCTGGCCCCCAGGCCCAGGTCAAGCGCCTTGGTTTGCCCACTAGGATTGTTTTAAGAAAATGGCAGACAAACCAGACATGGGGGAAATCGCCAGCTTCGATAAGGCCAAGCTGAAGAAAACGGAGACGCAGGAGAAGAACACCCTGCCGACCAAAGAGAGTGAGTGTGCCTCGGTCTCCCGCGCCCCAGCCCAGCCCCTCACCCTGCTCTTCCTTGCAAACCCACTCCTCCACCCCCCACCCCGCCGTTGTCCCCGGTGTGGGCGGCCCCGGCCACTCTTTCAGTTTCACAAAGCGCCTTGTTTCTCCCCAGCCCCAAGCTTCCTTCTAAATCCCCACACCTCGTGGGTGCCTCGCCCACACCGGGAAGCACCTCGGTTGCGGGTGGGGGTTGCAGCTCCCCTCCAGCGCCCGCTTCCCGCTCTCCACAGCCATTGAGCAGGAGAAGCGGAGTGAAATTTCCTAAGATCCTGGAGGATTTCCTACCCCCGTCCTCTTCGAGACCCCAGTCGTGATGTGGAGGAAGAGCCACCTGCAAGATGGACACGAGCCACAAGCTGCACTGTGAACCTGGGCACTCCGCGCCGATGCCACCGGCCTGTGGGTCTCTGAAGGGACCCCCCCCCAATCGGACTGCCAAATTCTCCGGTTTGCCCCGGGATATTATAGAAAATTATTTGTATGAATAATGAAAATAAAACACACCTCGTGGCA'

def test_custom_gene(dataframe_regression):
    cache = AssetCache()
    cache.set_custom_gene(name='USER_TMSB10', sequence=TMSB10_PREMRNA)

    data = generate_stub_data('USER_TMSB10', TMSB10_PREMRNA, first_n=10)

    df, aso_features = generate_aso_features(data, cache)

    dataframe_regression.check(df)

def test_short_gene(dataframe_regression, cache):
    target_gene = "TMSB10"  # Short gene
    transfection = Transfection.GYMNOSIS
    genome = "GRCh38"

    gene_to_data = cache.get_full_gene_data()
    gene_sequence = gene_to_data[target_gene].full_mrna

    data = generate_stub_data(
        target_gene=target_gene, gene_sequence=gene_sequence, first_n=10, transfection=transfection, genome=genome
    )
    df, aso_features = generate_aso_features(data, cache)

    # Check the resulting DataFrame against the saved regression file
    dataframe_regression.check(df)

#
# def test_short_gene_100(dataframe_regression, cache):
#     target_gene = "TMSB10"  # Short gene
#     transfection = Transfection.GYMNOSIS
#     genome = "GRCh38"
#
#     gene_to_data = cache.get_full_gene_data()
#     gene_sequence = gene_to_data[target_gene].full_mrna
#
#     data = generate_stub_data(
#         target_gene=target_gene, gene_sequence=gene_sequence, first_n=100, transfection=transfection, genome=genome
#     )
#     df, aso_features = generate_aso_features(data, cache)
#
#     # Check the resulting DataFrame against the saved regression file
#     dataframe_regression.check(df)
#
#
# def test_short_gene_1000(dataframe_regression, cache):
#     target_gene = "TMSB10"  # Short gene
#     genome = "GRCh38"
#     transfection = Transfection.GYMNOSIS
#
#     gene_to_data = cache.get_full_gene_data()
#     gene_sequence = gene_to_data[target_gene].full_mrna
#
#     data = generate_stub_data(
#         target_gene=target_gene, gene_sequence=gene_sequence, first_n=1000, transfection=transfection, genome=genome
#     )
#
#     df, aso_features = generate_aso_features(data, cache)
#
#     # Check the resulting DataFrame against the saved regression file
#     dataframe_regression.check(df)


def predict_aso_efficacy(model_path: str, input_df: pd.DataFrame) -> np.ndarray:
    """
    Loads a saved XGBoost model and predicts efficacy for the given DataFrame.
    Crashes (KeyError) if the DataFrame is missing required features.
    """
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")

    model = xgb.XGBRegressor()
    model.load_model(model_path)

    expected_cols = model.get_booster().feature_names
    if not expected_cols:
        raise ValueError(f"Model {model_path} does not contain embedded 'feature_names'.")

    # Strictly select expected columns. Will raise KeyError if any are missing.
    X = input_df[expected_cols].apply(pd.to_numeric, errors="coerce")

    # Predict (forces float32 for XGBoost compatibility)
    return model.predict(X.astype(np.float32).values)


def test_model_prediction(dataframe_regression, request):
    """
    Depends on the output of test_short_gene.
    Loads the saved CSV, runs inference, and checks for model regression.
    """
    # 1. Dynamically locate the CSV saved by test_short_gene
    # request.path gives the path to the current test file
    test_dir = request.path.parent
    test_filename_folder = request.path.stem
    csv_path = test_dir / test_filename_folder / "test_short_gene.csv"

    if not csv_path.exists():
        pytest.skip(f"Dependency missing: {csv_path} not found. Run test_short_gene first.")

    # 2. Load the features generated by the previous test
    df = pd.read_csv(csv_path)

    # 3. Path to your XGBoost JSON model
    # (Adjust this path if UnseenModel.json is stored elsewhere in your repo)
    model_path = NOTEBOOK_PATH / "models" / "UnseenOligoModel" / "UnseenModel.json"

    # 4. Run the prediction
    predictions = predict_aso_efficacy(model_path, df)

    # 5. Save the predictions as a DataFrame to track regression in model outputs
    predictions_df = pd.DataFrame({"Predicted_Inhibition": predictions})
    dataframe_regression.check(predictions_df)