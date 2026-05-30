"""Unit tests for cell-line name resolution (case- and punctuation-tolerant lookup)."""
import pytest

from tauso.data.consts import (
    CELL_LINE_TO_DEPMAP,
    CELL_LINE_TO_DEPMAP_PROXY_DICT,
    resolve_depmap_id,
    resolve_depmap_proxy,
    standardize_cell_line_name,
)


@pytest.mark.parametrize(
    "raw, expected_proxy",
    [
        # case variants of HeLa all collapse
        ("HeLa", "HeLa"),
        ("Hela", "HeLa"),
        ("HELA", "HeLa"),
        ("hela", "HeLa"),
        (" HeLa ", "HeLa"),
        # punctuation variants
        ("A431", "A-431"),
        ("A-431", "A-431"),
        ("HK-2", "HK-2"),
        ("HK2", "HK-2"),
        ("T-24", "T24"),
        ("T24", "T24"),
        ("SH-SY5Y", "SH-SY5Y"),
        ("SH-SY-5Y", "SH-SY5Y"),  # extra hyphen typo
        # period in dataset spelling
        ("MM.1R", "MM.1S"),
        ("MM.1S", "MM.1S"),
        ("mm1s", "MM.1S"),
    ],
)
def test_resolves_known_spelling_variants(raw, expected_proxy):
    assert resolve_depmap_proxy(raw) == expected_proxy


@pytest.mark.parametrize(
    "raw",
    [
        "HEK293T",  # explicitly blocked: too different from HEK293
        "HepaRG",  # primary-cell-derived; no honest cancer-line proxy
        "HUVEC",
        "iCell GABANeurons",
        "ReproNeuro Neurons (ReproCELL)",
        "KARPAS-229",  # B-cell line, not in DepMap; do not silently use Karpas-299 (T-cell ALCL)
    ],
)
def test_blocked_entries_resolve_to_none(raw):
    assert resolve_depmap_proxy(raw) is None
    assert resolve_depmap_id(raw) is None


@pytest.mark.parametrize("bad_input", [None, "", float("nan"), "zzz_unknown_line", 42])
def test_unknown_input_returns_none(bad_input):
    assert resolve_depmap_proxy(bad_input) is None
    assert resolve_depmap_id(bad_input) is None


@pytest.mark.parametrize(
    "raw, expected_id",
    [
        ("HeLa", "ACH-001086"),
        ("hela", "ACH-001086"),
        ("A431", "ACH-001328"),
        ("HepG2", "ACH-000739"),
        ("Hep3B", "ACH-000625"),
        ("U251", "ACH-000232"),  # U251 resolves through U-251 MG proxy
        ("MM.1R", "ACH-000763"),  # parental-line proxy
    ],
)
def test_resolve_depmap_id_composition(raw, expected_id):
    assert resolve_depmap_id(raw) == expected_id


@pytest.mark.parametrize(
    "raw, expected_tag",
    [
        # mapped: returns the proxy's stripped-uppercase tag
        ("HeLa", "HELA"),
        ("Hela", "HELA"),
        ("A-431", "A431"),
        # blocked: returns the cross-cell-line "Generic" fallback
        ("HepaRG", "Generic"),
        ("HUVEC", "Generic"),
        ("HEK293T", "Generic"),
        ("KARPAS-229", "Generic"),
        ("ReproNeuro Neurons (ReproCELL)", "Generic"),
        # empty / sentinel input
        ("None", "Generic"),
        ("none", "Generic"),
        (None, "Generic"),
        # unknown: fall back to stripped uppercase of raw input
        ("zzz_unknown_line", "ZZZUNKNOWNLINE"),
    ],
)
def test_standardize_cell_line_name(raw, expected_tag):
    assert standardize_cell_line_name(raw) == expected_tag


def test_proxy_dict_values_are_all_known_depmap_ids():
    """Every non-None proxy in CELL_LINE_TO_DEPMAP_PROXY_DICT must map to a real DepMap ID."""
    for input_name, proxy in CELL_LINE_TO_DEPMAP_PROXY_DICT.items():
        if proxy is None:
            continue
        assert proxy in CELL_LINE_TO_DEPMAP, f"Proxy {proxy!r} (from {input_name!r}) missing CELL_LINE_TO_DEPMAP entry"
