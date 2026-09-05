"""The {1st, 2nd, 3rd}-generation triple is read off MODIFICATION_STRING.

A gapmer's string names its high-affinity sugar, so the presence of "MOE" or "LNA"/"cEt" decides
the generation and their joint absence means a pure PS-DNA oligo. A mixmer's string is the bare
word "mixmer" and names no sugar at all, so that inference does not hold for it and the triple is
left unset.
"""

import numpy as np
import pandas as pd
import pytest

from tauso.data.consts import MIXMER_MODIFICATION, MODIFICATION_STRING, TRANSFECTION_RAW
from tauso.populate.calculators.calculator import Calculator

GENERATION_COLUMNS = ["chem_1st_gen", "chem_2nd_gen", "chem_3rd_gen"]


def _generations(modification_strings):
    data = pd.DataFrame(
        {
            "index_oligo": np.arange(len(modification_strings)),
            MODIFICATION_STRING: modification_strings,
            TRANSFECTION_RAW: ["Lipofection"] * len(modification_strings),
        }
    )
    calculator = Calculator(data=data)
    calculator.calculate_basic()
    return calculator.data[GENERATION_COLUMNS]


@pytest.mark.parametrize(
    "modification, expected",
    [
        ("DNA", [1, 0, 0]),
        ("MOE/5-methylcytosines/deoxy", [0, 1, 0]),
        ("cEt/5-methylcytosines/deoxy", [0, 0, 1]),
        ("LNA/5-methylcytosines/deoxy", [0, 0, 1]),
    ],
)
def test_gapmer_generations(modification, expected):
    assert _generations([modification]).iloc[0].tolist() == expected


def test_mixmer_leaves_every_generation_unset():
    """"mixmer" names no sugar, so it must not be read as the absence of one."""
    assert _generations([MIXMER_MODIFICATION]).iloc[0].isna().all()


def test_mixmer_does_not_disturb_its_neighbours():
    generations = _generations(["MOE/5-methylcytosines/deoxy", MIXMER_MODIFICATION, "DNA"])
    assert generations.iloc[0].tolist() == [0, 1, 0]
    assert generations.iloc[1].isna().all()
    assert generations.iloc[2].tolist() == [1, 0, 0]
