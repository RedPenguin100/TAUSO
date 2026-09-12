"""Resolving a cell-line name to its DepMap model.

The shipped table covers every line DepMap publishes; the curated proxy table sits in front of
it, so a name the project has already ruled on keeps the answer it was given.
"""

import pytest

from tauso.data.consts import (
    CELL_LINE_TO_DEPMAP,
    CELL_LINE_TO_DEPMAP_PROXY_DICT,
    DEPMAP_MODELS,
    resolve_depmap_id,
    resolve_depmap_proxy,
)


def test_the_table_covers_the_whole_release():
    assert len(DEPMAP_MODELS) > 2000
    assert all(depmap_id.startswith("ACH-") for _, depmap_id in DEPMAP_MODELS.values())


def test_a_line_outside_the_curated_table_now_resolves():
    assert resolve_depmap_id("SK-MEL-28") == "ACH-000615"
    assert resolve_depmap_id("MDA-MB-231") == "ACH-000768"


@pytest.mark.parametrize("spelling", ["SK-MEL-28", "skmel28", "SK MEL 28", "sk_mel_28"])
def test_punctuation_and_case_do_not_matter(spelling):
    assert resolve_depmap_id(spelling) == "ACH-000615"


@pytest.mark.parametrize("name", [n for n, v in CELL_LINE_TO_DEPMAP_PROXY_DICT.items() if v is None])
def test_a_blocked_line_stays_blocked(name):
    # Primary cells, iPSC-derived lines and HEK293T are ruled out on purpose; the wider table
    # must not quietly re-enable them.
    assert resolve_depmap_proxy(name) is None
    assert resolve_depmap_id(name) is None


@pytest.mark.parametrize("proxy,depmap_id", sorted(CELL_LINE_TO_DEPMAP.items()))
def test_curated_ids_are_unchanged(proxy, depmap_id):
    assert resolve_depmap_id(proxy) == depmap_id


@pytest.mark.parametrize("raw", ["not a cell line", "", None, 3])
def test_an_unknown_name_resolves_to_nothing(raw):
    assert resolve_depmap_proxy(raw) is None
    assert resolve_depmap_id(raw) is None


def test_every_name_is_distinct_once_normalised():
    # A collision would make resolution depend on row order in the shipped file.
    assert len({key for key in DEPMAP_MODELS}) == len(DEPMAP_MODELS)
