import pandas as pd
import pytest

from tauso.features.context.mrna_halflife import HALFLIFE_SPECIES, select_species

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
