"""
Verify that the _RC_RNA_TABLE str.translate approach produces byte-identical
FASTA content to the Bio.Seq.reverse_complement_rna() approach used in the
original get_trigger_mfe_scores_by_risearch.

If this test fails, the batch FASTA writer will produce wrong query sequences
and RIsearch will return empty/wrong results.
"""

import pytest
from Bio.Seq import Seq

from tauso.features.hybridization.fast_hybridization import _RC_RNA_TABLE
from tauso.util import get_antisense

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


def biopython_rc_rna(trigger: str) -> str:
    return str(Seq(trigger).reverse_complement_rna())


def translate_rc_rna(trigger: str) -> str:
    return trigger[::-1].translate(_RC_RNA_TABLE)


@pytest.mark.parametrize("aso_seq", SAMPLE_ASO_SEQS)
def test_rc_rna_equivalence(aso_seq):
    trigger = get_antisense(aso_seq)
    bio = biopython_rc_rna(trigger)
    fast = translate_rc_rna(trigger)
    assert fast == bio, f"Mismatch for aso_seq={aso_seq!r}: Bio.Seq gave {bio!r}, translate gave {fast!r}"


@pytest.mark.parametrize("aso_seq", SAMPLE_ASO_SEQS)
def test_fasta_block_content_matches(aso_seq):
    """The two-line FASTA block produced by each approach must be identical."""
    trigger = get_antisense(aso_seq)
    query_id = "42"

    bio_block = f">{query_id}\n{biopython_rc_rna(trigger)}\n"
    fast_block = f">{query_id}\n{translate_rc_rna(trigger)}\n"

    assert fast_block == bio_block


def test_multi_query_fasta_format():
    """Full multi-query FASTA built with the bulk-join approach equals line-by-line f.write."""
    pairs = [(str(i), get_antisense(seq)) for i, seq in enumerate(SAMPLE_ASO_SEQS)]

    # Line-by-line (original approach)
    old_parts = []
    for query_id, trigger in pairs:
        old_parts.append(f">{query_id}\n{biopython_rc_rna(trigger)}\n")
    old_content = "".join(old_parts)

    # Bulk join (new approach)
    lines = []
    for query_id, trigger in pairs:
        lines.append(f">{query_id}")
        lines.append(translate_rc_rna(trigger))
    new_content = "\n".join(lines) + "\n"

    assert new_content == old_content
