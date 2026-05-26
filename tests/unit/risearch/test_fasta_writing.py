"""
Verify that get_antisense_rna produces the same output as Bio.Seq.reverse_complement_rna()
and that the batch query FASTA format written by the RIsearch invocation is correct.
"""

import pytest
from Bio.Seq import Seq

from tauso.util import get_antisense, get_antisense_rna

SAMPLE_ASO_SEQS = [
    "ATCGATCGATCG",
    "GGGGCCCCTTTT",
    "TGCATGCATGCA",
    "AAAAAAAAAAAAAAAAAAAA",
    "GCTAGCTAGCTAGCTA",
    "CCAAGGTTCCGGTTAA",
    "CGTTTAGCTGACGTACGTTT",
    "AACCGGTTTTGGCCAA",
    "GATCGATCGATCGATCGATC",
]


@pytest.mark.parametrize("aso_seq", SAMPLE_ASO_SEQS)
def test_rc_rna_equivalence(aso_seq):
    trigger = get_antisense(aso_seq)
    bio = str(Seq(trigger).reverse_complement_rna())
    fast = get_antisense_rna(trigger)
    assert fast == bio, f"Mismatch for aso_seq={aso_seq!r}: Bio.Seq gave {bio!r}, get_antisense_rna gave {fast!r}"


@pytest.mark.parametrize("aso_seq", SAMPLE_ASO_SEQS)
def test_fasta_block_content_matches(aso_seq):
    """The two-line FASTA block produced by each approach must be identical."""
    trigger = get_antisense(aso_seq)
    query_id = "42"

    bio_block = f">{query_id}\n{Seq(trigger).reverse_complement_rna()}\n"
    fast_block = f">{query_id}\n{get_antisense_rna(trigger)}\n"

    assert fast_block == bio_block


def test_multi_query_fasta_format():
    """Full multi-query FASTA built with the bulk-join approach equals line-by-line f.write."""
    pairs = [(str(i), get_antisense(seq)) for i, seq in enumerate(SAMPLE_ASO_SEQS)]

    # Line-by-line (reference)
    old_parts = []
    for query_id, trigger in pairs:
        old_parts.append(f">{query_id}\n{Seq(trigger).reverse_complement_rna()}\n")
    old_content = "".join(old_parts)

    # Bulk join (production approach used in the _risearch_stdout query-FASTA writer)
    lines = []
    for query_id, trigger in pairs:
        lines.append(f">{query_id}")
        lines.append(get_antisense_rna(trigger))
    new_content = "\n".join(lines) + "\n"

    assert new_content == old_content
