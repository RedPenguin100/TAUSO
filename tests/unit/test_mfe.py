import numpy as np
import pandas as pd
import pytest

from tauso.features.fold.vienna_fold import calculate_avg_mfe_over_sense_region

# ---------------------------------------------------------
# Config
# ---------------------------------------------------------

SENSE = "ATGCGTACGTAGCTAGCTAG"
SENSE_LENGTH = len(SENSE)

FLANK_UNIT = "GCTA"

SETTINGS = [
    (120, 45, 7), (120, 45, 15), (120, 45, 10), (120, 45, 4),
    (120, 35, 15), (120, 35, 10), (120, 35, 7), (120, 35, 4),
    (120, 55, 15), (120, 55, 10), (120, 55, 7), (120, 55, 4),
    (120, 65, 15), (120, 65, 10), (120, 65, 7), (120, 65, 4),
]


def build_sequence(flank_size):
    flank = (FLANK_UNIT * ((flank_size // len(FLANK_UNIT)) + 1))[:flank_size]
    return flank + SENSE + flank


def test_mfe_regression_snapshot(dataframe_regression):

    rows = []

    for flank_size, window_size, step in SETTINGS:

        seq = build_sequence(flank_size)
        sense_start = flank_size

        for i in range(10):
            val = calculate_avg_mfe_over_sense_region(
                sequence=seq,
                sense_start_in_flank=sense_start,
                sense_length=SENSE_LENGTH,
                window_size=window_size,
                step=step,
            )

        rows.append({
            "flank_size": flank_size,
            "window_size": window_size,
            "step": step,
            "sequence_length": len(seq),
            "mfe": np.round(val, 6) if not np.isnan(val) else np.nan,
        })

    df = pd.DataFrame(rows).sort_values(
        ["flank_size", "window_size", "step"]
    ).reset_index(drop=True)

    dataframe_regression.check(df)