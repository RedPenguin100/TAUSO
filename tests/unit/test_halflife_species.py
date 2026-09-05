import pandas as pd
import pytest

from tauso.features.context.mrna_halflife import (
    HALFLIFE_BASELINE_CONDITIONS,
    HALFLIFE_SPECIES,
    select_conditions,
    select_species,
)

TTDB = pd.DataFrame(
    {
        "species_name": ["Human", "Mouse", " Human ", "Yeast"],
        "gene": ["ACTB", "Actb", "GAPDH", "ACT1"],
    }
)


def test_default_species_is_human():
    assert HALFLIFE_SPECIES == "Human"


def test_keeps_only_the_requested_species():
    assert select_species(TTDB, "Human").gene.tolist() == ["ACTB", "GAPDH"]
    assert select_species(TTDB, "Mouse").gene.tolist() == ["Actb"]


@pytest.mark.parametrize("spelling", ["human", "HUMAN", " Human "])
def test_species_name_matching_ignores_case_and_padding(spelling):
    assert select_species(TTDB, spelling).gene.tolist() == ["ACTB", "GAPDH"]


def test_absent_species_raises_and_names_what_is_available():
    with pytest.raises(ValueError) as e:
        select_species(TTDB, "Humann")
    assert "Human" in str(e.value) and "Mouse" in str(e.value)


CONDITIONS = pd.DataFrame(
    {
        "condition": ["WT", "uninfected", "heat shock time course", " FP control ", "DDX3X depletion"],
        "gene": ["ACTB", "GAPDH", "MYC", "TP53", "EGFR"],
    }
)


def test_default_conditions_are_the_baseline_arms():
    assert HALFLIFE_BASELINE_CONDITIONS[0] == "WT"
    assert "uninfected" in HALFLIFE_BASELINE_CONDITIONS


def test_keeps_baseline_arms_and_drops_treatments():
    kept = select_conditions(CONDITIONS, HALFLIFE_BASELINE_CONDITIONS)
    assert kept.gene.tolist() == ["ACTB", "GAPDH", "TP53"]


def test_absent_condition_raises_and_names_what_is_available():
    with pytest.raises(ValueError) as e:
        select_conditions(CONDITIONS, ["not a condition"])
    assert "heat shock time course" in str(e.value)
