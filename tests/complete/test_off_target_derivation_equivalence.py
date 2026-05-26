"""
Validate the two assumptions that let off-target derive several cutoffs / top_n subsets
from a single loose RIsearch run instead of one run each:

  1. cutoff: a stricter ``-s`` run equals the loose run filtered to ``score > cutoff``
     (RIsearch ``-s`` is exclusive, and the filter is on the score column, not energy).
  2. gene subset: running against all target genes and restricting to a gene subset
     equals running against just that subset (RIsearch scores each (query, target) pair
     independently).

Real genes, the real RIsearch binary, no mocks.
"""

import os

import pytest

from tauso.features.hybridization.fast_hybridization import (
    Interaction,
    dump_target_file,
    risearch_hits_dataframe,
)
from tauso.util import get_antisense

LOOSE, STRICT = 800, 1200
COLS = ["trigger", "target", "trigger_start", "target_start", "score", "energy"]
_RISEARCH = dict(parsing_type="2", interaction_type=Interaction.RNA_DNA_NO_WOBBLE, transpose=True)


@pytest.fixture(scope="module")
def target_and_queries(gene_to_data_full):
    """Three relatively short real genes as the off-target target, with one antisense
    query per gene. The query is a 120-mer so every gene's self-hit clears the strict
    cutoff regardless of its GC content; genes are kept >= 4000 nt so the region is well
    inside the transcript."""
    genes = sorted(
        (g for g, d in gene_to_data_full.items() if d.full_mrna and len(d.full_mrna) >= 4000),
        key=lambda g: len(gene_to_data_full[g].full_mrna),
    )[:3]
    seq_map = {g: gene_to_data_full[g].full_mrna for g in genes}
    queries = [(str(i), get_antisense(seq_map[g][500:620])) for i, g in enumerate(genes)]
    return genes, seq_map, queries


def _hits(df):
    return set(df[COLS].itertuples(index=False, name=None))


def test_cutoff_filter_uses_score_and_is_exclusive(target_and_queries):
    _, seq_map, queries = target_and_queries
    target = dump_target_file("ot-cutoff.fa", seq_map)
    try:
        loose = risearch_hits_dataframe(queries, target, minimum_score=LOOSE, **_RISEARCH)
        strict = risearch_hits_dataframe(queries, target, minimum_score=STRICT, **_RISEARCH)
    finally:
        os.remove(target)

    strict_hits = _hits(strict)
    assert strict_hits, "strict run produced no hits"
    assert strict_hits <= _hits(loose)  # loose is a superset of strict
    assert _hits(loose[loose["score"] > STRICT]) == strict_hits


def test_gene_subset_filter_matches_subset_target(target_and_queries):
    genes, seq_map, queries = target_and_queries
    subset = set(genes[:2])
    full_target = dump_target_file("ot-full.fa", seq_map)
    sub_target = dump_target_file("ot-sub.fa", {g: seq_map[g] for g in subset})
    try:
        full = risearch_hits_dataframe(queries, full_target, minimum_score=LOOSE, **_RISEARCH)
        sub = risearch_hits_dataframe(queries, sub_target, minimum_score=LOOSE, **_RISEARCH)
    finally:
        os.remove(full_target)
        os.remove(sub_target)

    sub_hits = _hits(sub)
    assert sub_hits, "subset run produced no hits"
    assert _hits(full[full["target"].isin(subset)]) == sub_hits


def test_cutoff_and_subset_filters_combine(target_and_queries):
    genes, seq_map, queries = target_and_queries
    subset = set(genes[:2])
    loose_full = dump_target_file("ot-comb-full.fa", seq_map)
    strict_sub = dump_target_file("ot-comb-sub.fa", {g: seq_map[g] for g in subset})
    try:
        loose = risearch_hits_dataframe(queries, loose_full, minimum_score=LOOSE, **_RISEARCH)
        direct = risearch_hits_dataframe(queries, strict_sub, minimum_score=STRICT, **_RISEARCH)
    finally:
        os.remove(loose_full)
        os.remove(strict_sub)

    direct_hits = _hits(direct)
    assert direct_hits, "direct run produced no hits"
    assert _hits(loose[(loose["score"] > STRICT) & (loose["target"].isin(subset))]) == direct_hits
