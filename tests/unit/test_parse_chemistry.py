"""The sugar pattern is positional, so it must have exactly one character per residue.

A residue the parser does not recognise used to contribute nothing, which shifted every
residue after it and left the pattern one short — silently, and with a label naming only the
gapmer sugar. These tests pin the length invariant and the mixmer labelling.
"""

import logging

import pytest
from notebooks.data.OligoAI.parse_chemistry import GAPMER_LABEL, SUGAR_CHAR, _process_chemistry

from tauso.data.consts import MIXMER_MODIFICATION

MOE_GAPMER = ["MOE"] * 5 + ["DNA"] * 10 + ["MOE"] * 5
CET_GAPMER = ["CET"] * 3 + ["DNA"] * 10 + ["CET"] * 3


def _parse(mods):
    return _process_chemistry(str(mods))


@pytest.mark.parametrize(
    "mods, pattern, label",
    [
        (MOE_GAPMER, "MMMMMddddddddddMMMMM", "MOE/5-methylcytosines/deoxy"),
        (CET_GAPMER, "CCCddddddddddCCC", "cEt/5-methylcytosines/deoxy"),
        (["LNA", "DNA", "LNA"], "LdL", "LNA/5-methylcytosines/deoxy"),
        (["DNA"] * 4, "dddd", "DNA"),
        # more than one gapmer sugar
        (["MOE", "CET", "DNA"], "MCd", "mixmer"),
    ],
)
def test_supported_chemistries_are_unchanged(mods, pattern, label):
    assert _parse(mods) == (pattern, label)


@pytest.mark.parametrize("other", ["OME", "F"])
def test_a_sugar_tauso_does_not_model_makes_it_a_mixmer(other):
    """2'-OMe and 2'-F are real chemistries, not gapmer sugars. They must not be mistaken for
    a plain cEt gapmer — the label is what excludes the row upstream."""
    mods = ["CET", "CET", "CET", "DNA", other] + ["DNA"] * 9 + ["CET", "CET", "CET"]
    pattern, label = _parse(mods)

    assert label == "mixmer"
    assert len(pattern) == len(mods), "one character per residue, or positions shift"


@pytest.mark.parametrize("mods", [MOE_GAPMER, CET_GAPMER, ["CET", "OME", "DNA"], ["F", "DNA"]])
def test_pattern_has_one_character_per_residue(mods):
    pattern, _ = _parse(mods)
    assert len(pattern) == len(mods)


def test_an_unrecognised_sugar_is_refused_not_dropped(caplog):
    """The failure that started this: an unmapped residue must never yield a short pattern."""
    mods = ["CET", "GLYCOL", "DNA"]
    with caplog.at_level(logging.WARNING):
        pattern, label = _parse(mods)

    assert (pattern, label) == (None, None)
    assert "GLYCOL" in caplog.text


def test_every_mapped_sugar_has_a_distinct_character():
    assert len(set(SUGAR_CHAR.values())) == len(SUGAR_CHAR)


def test_unparsable_input_is_still_tolerated():
    assert _process_chemistry(None) == (None, None)
    assert _process_chemistry("not a list") == (None, None)
    assert _process_chemistry("[unclosed") == (None, None)


def test_the_mixmer_label_is_what_excludes_the_row():
    """The parser only labels; preprocessing does the dropping. Pin the link between them so a
    rename on either side cannot quietly let an unmodelled sugar through."""
    from notebooks.preprocessing import SUPPORTED_CHEMISTRIES

    assert MIXMER_MODIFICATION not in SUPPORTED_CHEMISTRIES
    for label in GAPMER_LABEL.values():
        assert label in SUPPORTED_CHEMISTRIES
