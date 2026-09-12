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

import pytest
from pyrisearch_tauso import fasta_targets

from tauso.genome.read_human_genome import get_locus_to_data_dict

from .hits import risearch_hits_dataframe

LOOSE, STRICT = 800, 1200
COLS = ["query", "target", "query_start", "target_start", "score", "energy"]
_RISEARCH = dict(transpose=True)


@pytest.fixture(scope="module")
def gene_to_data_full():
    """Locus-to-data dict for all genes, with introns."""
    return get_locus_to_data_dict(include_introns=True)


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
    queries = [(str(i), seq_map[g][500:620]) for i, g in enumerate(genes)]
    return genes, seq_map, queries


def _hits(df):
    return set(df[COLS].itertuples(index=False, name=None))


def test_cutoff_filter_uses_score_and_is_exclusive(target_and_queries):
    _, seq_map, queries = target_and_queries
    with fasta_targets(seq_map) as target:
        loose = risearch_hits_dataframe(queries, target, minimum_score=LOOSE, **_RISEARCH)
        strict = risearch_hits_dataframe(queries, target, minimum_score=STRICT, **_RISEARCH)

    strict_hits = _hits(strict)
    assert strict_hits, "strict run produced no hits"
    assert strict_hits <= _hits(loose)  # loose is a superset of strict
    assert _hits(loose[loose["score"] > STRICT]) == strict_hits


def test_gene_subset_filter_matches_subset_target(target_and_queries):
    genes, seq_map, queries = target_and_queries
    subset = set(genes[:2])
    with fasta_targets(seq_map) as full_target, fasta_targets({g: seq_map[g] for g in subset}) as sub_target:
        full = risearch_hits_dataframe(queries, full_target, minimum_score=LOOSE, **_RISEARCH)
        sub = risearch_hits_dataframe(queries, sub_target, minimum_score=LOOSE, **_RISEARCH)

    sub_hits = _hits(sub)
    assert sub_hits, "subset run produced no hits"
    assert _hits(full[full["target"].isin(subset)]) == sub_hits


def test_cutoff_and_subset_filters_combine(target_and_queries):
    genes, seq_map, queries = target_and_queries
    subset = set(genes[:2])
    with fasta_targets(seq_map) as loose_full, fasta_targets({g: seq_map[g] for g in subset}) as strict_sub:
        loose = risearch_hits_dataframe(queries, loose_full, minimum_score=LOOSE, **_RISEARCH)
        direct = risearch_hits_dataframe(queries, strict_sub, minimum_score=STRICT, **_RISEARCH)

    direct_hits = _hits(direct)
    assert direct_hits, "direct run produced no hits"
    assert _hits(loose[(loose["score"] > STRICT) & (loose["target"].isin(subset))]) == direct_hits
