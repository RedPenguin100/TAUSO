"""A second `set_custom_gene` for the same gene name must not keep the first sequence.

The lean cache and the gene registry are keyed on the gene NAMES, which do not change when
only the sequence does, so nothing in the key tells them the data went stale.
"""

import pytest

import tauso.genome.TranscriptMapper as mapper
from tauso.populate.calculators import cache as cache_module
from tauso.populate.calculators.cache import AssetCache

GENE = "CUSTOM"
FIRST = "ACGT" * 10
SECOND = "TTTT" * 10


@pytest.fixture
def asset_cache(monkeypatch):
    """An AssetCache whose genome lookup returns nothing, so only the caching logic runs."""
    monkeypatch.setattr(cache_module, "get_locus_to_data_dict", lambda **kwargs: {})
    return AssetCache()


def test_second_custom_gene_replaces_the_first_sequence(asset_cache):
    asset_cache.set_custom_gene(GENE, FIRST)
    assert asset_cache.get_lean_gene({GENE})[GENE].full_mrna == FIRST

    asset_cache.set_custom_gene(GENE, SECOND)
    assert asset_cache.get_lean_gene({GENE})[GENE].full_mrna == SECOND


def test_the_registry_does_not_serve_the_first_sequence_either(asset_cache, monkeypatch):
    seen = []
    monkeypatch.setattr(
        mapper, "build_gene_sequence_registry", lambda genes, gene_to_data: seen.append(gene_to_data[GENE].full_mrna)
    )
    asset_cache.set_custom_gene(GENE, FIRST)
    asset_cache.get_gene_registry({GENE})
    asset_cache.set_custom_gene(GENE, SECOND)
    asset_cache.get_gene_registry({GENE})

    assert seen == [FIRST, SECOND]


def test_cds_annotation_is_replaced_too(asset_cache):
    asset_cache.set_custom_gene(GENE, FIRST, cds_start=0, cds_end=12)
    asset_cache.get_lean_gene({GENE})  # populate the cache, so the second call has one to drop
    asset_cache.set_custom_gene(GENE, SECOND, cds_start=3, cds_end=21)
    locus = asset_cache.get_lean_gene({GENE})[GENE]
    assert locus.full_mrna == SECOND
    assert (locus.start_codon, locus.stop_codon) == ((3, 6), (18, 21))
    assert locus.five_prime_utr == SECOND[:3]
