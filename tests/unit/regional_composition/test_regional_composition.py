"""Dinucleotide composition per gapmer region, measured against the rest of the oligo.

These pin the structure rather than the numbers: the regions partition the steps exactly, a
comparison that cannot be made is missing rather than zero, and every region's values cancel.
"""

import numpy as np
import pandas as pd
import pytest

from tauso.features.regional_composition.regional_composition import (
    DINUCLEOTIDES,
    FEATURE_NAMES,
    REGIONS,
    calculate_regional_composition,
    step_regions,
)
from tauso.populate.populate_regional_composition import (
    populate_regional_composition_features,
    regional_composition_feature_name,
    regional_composition_feature_names,
)

MOE_GAPMER = "MMMMMddddddddddMMMMM"
CET_GAPMER = "CCCddddddddddCCC"
ALL_DNA = "dddddddddddddddd"
SEQ_16 = "AAACCTAAAATAGTGG"
SEQ_20 = "AAAAAAAAACCTAATAGACG"


def score(sequence, pattern):
    return {k: float(v[0]) for k, v in calculate_regional_composition([sequence], [pattern]).items()}


# --- the regions partition the steps ------------------------------------------------
@pytest.mark.parametrize(
    "pattern,length",
    [
        (MOE_GAPMER, 20),
        (CET_GAPMER, 16),
        (ALL_DNA, 16),
        ("MMMMddddddddddMMMMMd", 20),
        ("dCddCddCdddddddd", 16),
        ("ddddddddddMMMMMM", 16),
        ("MMMMMMdddddddddd", 16),
    ],
)
def test_regions_cover_every_step_exactly_once(pattern, length):
    masks = np.vstack([step_regions(pattern, length)[r] for r in REGIONS])
    assert masks.sum() == length - 1
    assert masks.sum(axis=0).max() <= 1


def test_wing_grows_with_the_wing_length():
    short = step_regions(CET_GAPMER, 16)
    long = step_regions(MOE_GAPMER, 20)
    assert short["wing5"].sum() == 2
    assert long["wing5"].sum() == 4
    assert short["gap"].sum() == long["gap"].sum() == 9


def test_junction_is_a_single_step_on_a_plain_gapmer():
    masks = step_regions(MOE_GAPMER, 20)
    assert masks["j5"].sum() == 1
    assert masks["j3"].sum() == 1


# --- what the values must satisfy ---------------------------------------------------
def test_values_cancel_within_each_region():
    # Both the region and its reference are frequencies over sixteen cells, so they sum to one.
    scored = score(SEQ_20, MOE_GAPMER)
    for region in REGIONS:
        total = sum(scored[f"{d}_{region}"] for d in DINUCLEOTIDES)
        assert total == pytest.approx(0.0, abs=1e-12)


def test_an_all_dna_oligo_has_nothing_to_compare():
    # The gap holds every step and the wings hold none, so no region has a reference.
    scored = score(SEQ_16, ALL_DNA)
    assert all(np.isnan(v) for v in scored.values())


def test_an_empty_region_is_missing_not_zero():
    scored = score(SEQ_16, "MMMMMMMMMMMMdddd")
    assert all(np.isnan(scored[f"{d}_wing3"]) for d in DINUCLEOTIDES)


def test_a_region_is_not_inside_its_own_reference():
    # AAACCTAAAATAGTGG on MMMddddddddddMMM has AA on 5 of its 15 steps, 2 of them the 5' wing.
    # The wing is read against the other 13 steps, which hold the remaining 3.
    scored = score(SEQ_16, "MMMddddddddddMMM")
    assert scored["AA_wing5"] == pytest.approx(1.0 - 3 / 13)
    # AC and GT each occur once, at the two junction steps, so nothing is left outside.
    assert scored["AC_j5"] == pytest.approx(1.0)
    assert scored["GT_j3"] == pytest.approx(1.0)


def test_chemistry_moves_the_boundary():
    moved = score(SEQ_20, "MMMMMMMdddddddddddd")
    plain = score(SEQ_20, MOE_GAPMER)
    assert any(not np.isclose(moved[f"{d}_gap"], plain[f"{d}_gap"], equal_nan=True) for d in DINUCLEOTIDES)


# --- the populate step ---------------------------------------------------------------
def test_populate_adds_every_column():
    df = pd.DataFrame({"aso_sequence": [SEQ_16, SEQ_20], "chemical_pattern": [CET_GAPMER, MOE_GAPMER]})
    out, added = populate_regional_composition_features(df)
    assert added == regional_composition_feature_names()
    assert len(added) == 80
    assert all(c in out.columns for c in added)


def test_populate_needs_the_chemistry():
    with pytest.raises(ValueError):
        populate_regional_composition_features(pd.DataFrame({"aso_sequence": [SEQ_16]}))


def test_feature_name_rejects_an_unknown_quantity():
    with pytest.raises(ValueError):
        regional_composition_feature_name("ZZ_wing5")


def test_names_are_the_sixteen_dinucleotides_in_five_regions():
    assert len(FEATURE_NAMES) == len(set(FEATURE_NAMES)) == 80
    assert sorted(FEATURE_NAMES) == sorted(f"{d}_{r}" for r in REGIONS for d in DINUCLEOTIDES)
