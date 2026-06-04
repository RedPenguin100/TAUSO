"""CanonicalCds: coding-exon = exon - UTR, and correct handling on both strands.

The interval arithmetic is strand-agnostic (genomic coordinates); the strand is resolved
only when coding spans are converted to pre-mRNA coordinates. These tests pin both halves,
including a negative-strand end-to-end check (BRCA1 / ACTB / TP53 are minus-strand).
"""

import pytest

from tauso.genome.LocusInfo import CanonicalCds, StrandType, _coding_exon_intervals
from tauso.genome.read_human_genome import get_locus_to_data_dict

_START_CODONS = {"ATG", "AUG"}


# --- pure interval arithmetic (no genome needed, strand-agnostic) ---


def test_coding_exon_intervals_removes_terminal_utrs():
    # exon1 carries a 5' UTR head, exon2 a 3' UTR tail; the CDS is what is left.
    assert _coding_exon_intervals([(0, 50), (100, 160)], [(0, 20), (150, 160)]) == [(20, 50), (100, 150)]


def test_coding_exon_intervals_drops_fully_utr_exon():
    assert _coding_exon_intervals([(0, 30), (40, 60)], [(0, 30)]) == [(40, 60)]


def test_coding_exon_intervals_no_utr_is_identity():
    assert _coding_exon_intervals([(5, 10), (20, 25)], []) == [(5, 10), (20, 25)]


def test_coding_exon_intervals_internal_utr_splits_exon():
    assert _coding_exon_intervals([(0, 100)], [(40, 60)]) == [(0, 40), (60, 100)]


# --- end-to-end on real genes (both strands) ---


@pytest.fixture(scope="module")
def canonical():
    return get_locus_to_data_dict(gene_subset=["ACTB", "BRCA1", "TP53", "GAPDH"], canonical_only=True)


@pytest.mark.parametrize("gene", ["ACTB", "BRCA1", "TP53"])
def test_negative_strand_cds_valid_and_maps_start_codon(canonical, gene):
    info = canonical[gene]
    assert info.strand == StrandType.NEG
    cds = CanonicalCds.from_locus(info)
    # The biological 5' coding base maps to CDS index 0 (the start codon).
    first_premrna_start = cds.premrna_intervals[0][0]
    assert cds.to_cds_index(first_premrna_start) == 0
    assert cds.sequence[:3] in _START_CODONS
    # The pre-mRNA spans reconstruct the CDS in order (contiguous across splice joins).
    assert "".join(info.full_mrna[s:e] for s, e, _ in cds.premrna_intervals) == cds.sequence


def test_positive_strand_cds_maps_start_codon(canonical):
    info = canonical["GAPDH"]
    assert info.strand == StrandType.POS
    cds = CanonicalCds.from_locus(info)
    assert cds.to_cds_index(cds.premrna_intervals[0][0]) == 0
    assert cds.sequence[:3] in _START_CODONS
