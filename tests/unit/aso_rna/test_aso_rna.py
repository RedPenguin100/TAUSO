"""Geometry of the ASO:RNA duplex, averaged over the wings and the gap.

These pin the structure rather than the numbers: which table cell a step is scored in, that
the chemistry changes the answer, and that a region with nothing in it is missing. The numbers
themselves are pinned in tests/complete/regression.
"""

import numpy as np
import pandas as pd
import pytest

from tauso.features.aso_rna.aso_rna import (
    FEATURE_NAMES,
    JUNCTION_MEAN,
    OBSERVABLES,
    REGIONS,
    UNIFORM_MEAN,
    UNIFORM_SPREAD,
    calculate_aso_rna,
    single_dna_gap,
    step_cell,
    step_regions,
    sugars,
)
from tauso.populate.populate_aso_rna import (
    aso_rna_feature_name,
    aso_rna_feature_names,
    populate_aso_rna_features,
)

MOE_GAPMER = "MMMMMddddddddddMMMMM"
CET_GAPMER = "CCCddddddddddCCC"
ALL_DNA = "dddddddddddddddd"
SEQ_16 = "AAACCTAAAATAGTGG"
SEQ_20 = "AAAAAAAAACCTAATAGACG"


def score(sequence, pattern):
    return {k: float(v[0]) for k, v in calculate_aso_rna([sequence], [pattern]).items()}


# --- which cell a step is scored in --------------------------------------------------
def test_a_uniform_step_uses_the_sugar_against_rna():
    assert step_cell("M", "M", "AC") == ("M:R@AC", False)
    assert step_cell("D", "D", "AC") == ("D:R@AC", False)


def test_a_step_that_changes_chemistry_uses_the_junction_cell():
    assert step_cell("M", "D", "AC") == ("MD|RR@AC", True)
    assert step_cell("D", "E", "AC") == ("DE|RR@AC", True)


def test_cet_is_read_as_the_letter_the_tables_use():
    # The chemical_pattern column writes cEt as "C"; the weight tables name it "E".
    assert sugars(CET_GAPMER, 16)[0] == "E"
    assert step_cell("E", "D", "AA")[0] in JUNCTION_MEAN["Roll"]


def test_every_junction_cell_a_gapmer_needs_is_present():
    for sugar in ("M", "E"):
        for dinucleotide in ("AA", "CG", "TT"):
            assert step_cell(sugar, "D", dinucleotide)[0] in JUNCTION_MEAN["Roll"]
            assert step_cell("D", sugar, dinucleotide)[0] in JUNCTION_MEAN["Roll"]


@pytest.mark.parametrize("pattern", ["CCCdoddddddddCCC", "CCCdfddddddddCCC", "LLLddddddddddLLL", "CCCdxddddddddCCC"])
def test_a_sugar_the_tables_do_not_cover_makes_the_oligo_unscorable(pattern):
    # 2'-O-methyl, 2'-fluoro and LNA all appear in chemical_pattern; none has cells here, and
    # reading them as deoxy would be worse than declining to score.
    assert all(np.isnan(v) for v in score(SEQ_16, pattern).values())


@pytest.mark.parametrize(
    "pattern,why",
    [
        ("MMMddMddMddMddMddMMM", "five tied deoxy runs, so which is the gap turns on scan order"),
        ("MMMMddddddddddMMMMMd", "a stray deoxy in the 3' wing"),
        ("M" * 20, "no deoxy at all"),
    ],
)
def test_without_one_deoxy_stretch_there_is_nothing_to_measure(pattern, why):
    assert single_dna_gap(pattern) is None, why
    assert all(np.isnan(v) for v in score("A" * 20, pattern).values())


def test_one_deoxy_stretch_is_the_gap():
    assert single_dna_gap(MOE_GAPMER) == (5, 15)
    assert single_dna_gap(ALL_DNA) == (0, 16)


def test_a_pattern_that_does_not_match_the_sequence_raises():
    with pytest.raises(ValueError):
        calculate_aso_rna([SEQ_16], [CET_GAPMER[:-1]])


def test_the_three_covered_sugars_score():
    for pattern in ("CCCddddddddddCCC", "MMMddddddddddMMM", "d" * 16):
        assert any(np.isfinite(v) for v in score(SEQ_16, pattern).values())


def test_the_tables_carry_every_observable():
    for tables in (UNIFORM_MEAN, UNIFORM_SPREAD, JUNCTION_MEAN):
        assert set(tables) == set(OBSERVABLES)
        assert all(tables[o] for o in OBSERVABLES)


# --- the regions ----------------------------------------------------------------------
def test_a_boundary_step_belongs_to_no_region():
    masks = step_regions(MOE_GAPMER, 20)
    covered = sum(int(m.sum()) for m in masks.values())
    assert covered == 19 - 2


def test_the_wing_grows_with_the_wing_length():
    assert step_regions(CET_GAPMER, 16)["wing5"].sum() == 2
    assert step_regions(MOE_GAPMER, 20)["wing5"].sum() == 4


def test_an_oligo_with_no_gap_is_missing_everywhere():
    scored = score(SEQ_20, "M" * 20)
    assert all(np.isnan(v) for v in scored.values())


def test_an_all_dna_oligo_has_no_wings():
    scored = score(SEQ_16, ALL_DNA)
    assert all(np.isnan(scored[f"{o.lower()}_wing5"]) for o in OBSERVABLES)
    assert all(np.isnan(scored[f"{o.lower()}_wing3"]) for o in OBSERVABLES)
    assert all(np.isfinite(scored[f"{o.lower()}_gap"]) for o in OBSERVABLES)


# --- what the values must satisfy -----------------------------------------------------
def test_chemistry_changes_the_answer():
    moe = score(SEQ_16, "MMMddddddddddMMM")
    cet = score(SEQ_16, CET_GAPMER)
    assert any(not np.isclose(moe[f"{o.lower()}_wing5"], cet[f"{o.lower()}_wing5"]) for o in OBSERVABLES)


def test_the_gap_is_read_against_rna_as_dna():
    # Every gap step is deoxy on both sides, so the gap mean is a mean of D:R cells.
    scored = score(SEQ_20, MOE_GAPMER)
    cells = [UNIFORM_MEAN["Roll"][f"D:R@{SEQ_20[i - 1] + SEQ_20[i]}"] for i in range(6, 15)]
    assert scored["roll_gap"] == pytest.approx(float(np.mean(cells)))


def test_spread_and_mean_are_different_readings():
    scored = score(SEQ_20, MOE_GAPMER)
    assert any(not np.isclose(scored[f"{o.lower()}_gap"], scored[f"{o.lower()}_gap_spread"]) for o in OBSERVABLES)


# --- the populate step ----------------------------------------------------------------
def test_populate_adds_every_column():
    df = pd.DataFrame({"aso_sequence": [SEQ_16, SEQ_20], "chemical_pattern": [CET_GAPMER, MOE_GAPMER]})
    out, added = populate_aso_rna_features(df)
    assert added == aso_rna_feature_names()
    assert len(added) == 78
    assert all(c in out.columns for c in added)


def test_populate_needs_the_chemistry():
    with pytest.raises(ValueError):
        populate_aso_rna_features(pd.DataFrame({"aso_sequence": [SEQ_16]}))


def test_feature_name_rejects_an_unknown_quantity():
    with pytest.raises(ValueError):
        aso_rna_feature_name("roll_middle")


def test_names_are_every_observable_in_every_region_twice():
    assert len(FEATURE_NAMES) == len(set(FEATURE_NAMES)) == len(OBSERVABLES) * len(REGIONS) * 2
