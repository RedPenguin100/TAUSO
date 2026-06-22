"""Calculate features for the ASOptimizer dataset.

Loads ``PROCESSED_CSV``, standardizes the ASOptimizer columns to TAUSO's canonical
schema, builds a Calculator (``data_version='asoptimizer'``) and runs the full feature
pipeline (``Calculator.calculate_all``), writing the features to the feature store.

Usage:
    python -m notebooks.data.ASOptimizer.asoptimizer_features_calculate --cpus 32
"""
import argparse
import logging

import pandas as pd
from notebooks.data.ASOptimizer.consts import PROCESSED_CSV
from notebooks.data.ASOptimizer.utility import standardize_asoptimizer_data

from tauso.populate.calculators.calculator import Calculator

logger = logging.getLogger(__name__)

DATA_VERSION = "asoptimizer"


def load_data():
    logger.info("loading %s", PROCESSED_CSV)
    df = pd.read_csv(PROCESSED_CSV)
    logger.info("loaded %d rows", len(df))
    return standardize_asoptimizer_data(df)


def calculate_features(df, cpus=32, overwrite=True):
    calculator = Calculator(data=df, data_version=DATA_VERSION, overwrite=overwrite, cpus=cpus)
    calculator.calculate_all()
    return calculator


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpus", type=int, default=32, help="CPUs for parallel feature steps")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s | %(message)s",
    )

    df = load_data()
    calculate_features(df, cpus=args.cpus)
    logger.info("done")


if __name__ == "__main__":
    main()
