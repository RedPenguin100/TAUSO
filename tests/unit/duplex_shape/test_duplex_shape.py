"""Duplex geometry by region: what must hold whatever the weight table says."""

import numpy as np
import pandas as pd
import pytest

from tauso.features.duplex_shape.duplex_shape import (
    OBSERVABLES,
    QUANTITIES,
    REGIONS,
    TABLES,
    calculate_duplex_shape,
    step_masks,
    sugars,
)
from tauso.populate.populate_duplex_shape import (
    duplex_shape_feature_name,
    duplex_shape_feature_names,
    populate_duplex_shape_features,
)

MOE_GAPMER = "MMMMMddddddddddMMMMM"
CET_GAPMER = "CCCddddddddddddddCCC"
SEQ = "GCATTGGTACCTAGCAGTCA"


def score(sequence, pattern):
    return {k: float(v[0]) for k, v in calculate_duplex_shape([sequence], [pattern]).items()}


def test_one_quantity_per_observable_and_region():
    assert len(QUANTITIES) == len(OBSERVABLES) * len(REGIONS)
    assert set(QUANTITIES) == {f"{o}_{r}" for o in OBSERVABLES for r in REGIONS}


def test_every_table_covers_every_dinucleotide_in_every_state():
    for name, table in TABLES.items():
        states = {state for state, _ in table}
        assert states == {"D:D", "D:E", "D:M", "E:E", "M:M"}, name
        for state in states:
            cells = {cell for s, cell in table if s == state}
            assert len(cells) == 16, f"{name} {state}: {len(cells)} dinucleotides"


def test_regions_follow_the_dna_gap():
    masks = step_masks(MOE_GAPMER, len(MOE_GAPMER))
    steps = np.arange(1, len(MOE_GAPMER))
    # wings are positions 0-4 and 15-19, so the steps wholly inside them are 1-4 and 16-19
    assert set(steps[masks["wing5"]]) == {1, 2, 3, 4}
    assert set(steps[masks["wing3"]]) == {16, 17, 18, 19}
    assert set(steps[masks["gap"]]) == set(range(6, 15))


def test_boundary_steps_belong_to_no_region():
    """The two steps straddling a wing/gap boundary are counted nowhere."""
    masks = step_masks(MOE_GAPMER, len(MOE_GAPMER))
    for step in (5, 15):
        index = step - 1
        assert not any(mask[index] for mask in masks.values()), step


def test_regions_partition_without_overlap():
    masks = step_masks(MOE_GAPMER, len(MOE_GAPMER))
    total = sum(mask.astype(int) for mask in masks.values())
    assert total.max() <= 1


def test_short_wing_still_has_wing_steps():
    """A 3-nt wing leaves two steps inside it; nothing is dropped from the ends."""
    masks = step_masks("MMMddddddddddMMM", 16)
    assert masks["wing5"].sum() == 2
    assert masks["wing3"].sum() == 2


def test_wings_are_nan_without_a_wing():
    out = score("ACGT" * 5, "d" * 20)
    for observable in OBSERVABLES:
        assert np.isnan(out[f"{observable}_wing5"])
        assert np.isnan(out[f"{observable}_wing3"])
        assert not np.isnan(out[f"{observable}_gap"])


def test_chemistry_changes_the_wings_but_not_the_gap():
    moe = score(SEQ, MOE_GAPMER)
    cet = score(SEQ, "CCCCCddddddddddCCCCC")
    assert moe["roll_wing5"] != cet["roll_wing5"]
    assert moe["roll_gap"] == pytest.approx(cet["roll_gap"])


def test_asymmetric_wings_are_scored_separately():
    out = score(SEQ, "MMMMMMMMdddddddddMMM")
    assert not np.isnan(out["roll_wing5"])
    assert not np.isnan(out["roll_wing3"])
    assert out["roll_wing5"] != out["roll_wing3"]


def test_unknown_pattern_letter_is_deoxy():
    assert sugars("MMxdd", 5) == ["M", "M", "D", "D", "D"]


def test_cet_encoded_from_C():
    assert sugars(CET_GAPMER, len(CET_GAPMER))[0] == "E"


def test_populate_columns():
    df = pd.DataFrame({"aso_sequence": [SEQ, SEQ], "chemical_pattern": [MOE_GAPMER, CET_GAPMER]})
    out, columns = populate_duplex_shape_features(df)
    assert columns == duplex_shape_feature_names()
    assert all(c in out.columns for c in columns)
    assert len(out) == 2


def test_populate_missing_input():
    with pytest.raises(ValueError, match="Missing columns"):
        populate_duplex_shape_features(pd.DataFrame({"aso_sequence": [SEQ]}))


def test_unknown_quantity_rejected():
    with pytest.raises(ValueError, match="unknown quantity"):
        duplex_shape_feature_name("roll_middle")
