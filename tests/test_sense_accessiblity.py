import pytest
import numpy as np

import pandas as pd

from tauso.features.rna_access.rna_access import RNAAccess, parse
from tauso.new_model.consts_dataframe import *
from tauso.features.feature_names import SENSE_START, SENSE_LENGTH, SENSE_START_FROM_END
from tauso.new_model.populate.populate_sense_accessibility import populate_sense_accessibility, SENSE_AVG_ACCESSIBILITY

from tauso.new_model.consts_dataframe import SEQUENCE
from tauso.util import get_antisense

FLANK_SIZE = 120
ACCESS_SIZE = 13
SEED_SIZE = 13
SEED_SIZES = [SEED_SIZE * m for m in range(1, 4)]
ACCESS_WIN_SIZE = 80


def get_init_df(target_mrna, aso_sizes=[16, 20], canonical_name='DDX11L1'):
    candidates = []
    sense_starts = []
    sense_lengths = []
    sense_starts_from_end = []

    for aso_size in aso_sizes:
        for i in range(0, len(target_mrna) - (aso_size - 1)):
            target = target_mrna[i: i + aso_size]
            candidates.append(get_antisense(str(target)))
            sense_starts.append(i)
            sense_lengths.append(aso_size)
            sense_starts_from_end.append(i)
    df = pd.DataFrame({SEQUENCE: candidates, SENSE_START: sense_starts, SENSE_LENGTH: sense_lengths,
                       SENSE_START_FROM_END: sense_starts_from_end, CANONICAL_GENE: canonical_name})
    return df


def test_regression(short_mrna, n_compare=10, path="avg_sense_access.txt", use_saved=True):
    df = get_init_df(short_mrna.full_mrna, [16])
    print("Length mRNA: ", len(short_mrna.full_mrna))

    populate_sense_accessibility(df, short_mrna)
    avg_sense_predictions = list(df[SENSE_AVG_ACCESSIBILITY])

    if not use_saved:
        # write predictions so they can serve as regression data
        with open(path, "w") as f:
            for val in avg_sense_predictions:
                f.write(f"{val if val is not None else 'nan'}\n")

            avg_sense_regression = np.array(avg_sense_predictions)
    else:
        avg_sense_regression = np.loadtxt(path)
    # compare
    np.testing.assert_allclose(
        np.array(avg_sense_predictions[:n_compare]),
        np.array(avg_sense_regression[:n_compare])
    )



def test_parsing(dataframe_regression):
    with open("trig_raccess.txt", "r") as f:
        data = f.read()

    for i in range(10000):
        seed_sizes = [SEED_SIZE]
        res = parse(data, seed_sizes)

    # Expect a single DataFrame under the key 'rna'
    dataframe_regression.check(res["rna"])