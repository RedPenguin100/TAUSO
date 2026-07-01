"""Optional CDS annotation for custom genes (LocusInfo / set_custom_gene)."""

import pytest

from tauso.genome.LocusInfo import CanonicalCds, GeneType, LocusInfo

SEQ = "ACGTACGT" * 20  # 160 nt
CDS_START, CDS_END = 10, 40  # 30 nt in-frame coding span


def test_no_cds_is_unannotated():
    loc = LocusInfo(seq=SEQ)
    assert loc.gene_type == GeneType.UNKNOWN
    assert loc._5utr_indices == [] and loc._3utr_indices == []
    assert CanonicalCds.from_locus(loc).sequence == ""  # non-coding -> empty CDS


def test_cds_annotation():
    loc = LocusInfo(seq=SEQ, cds_start=CDS_START, cds_end=CDS_END)
    assert loc.gene_type == GeneType.PROTEIN_CODING
    assert loc._5utr_indices == [(0, CDS_START)]
    assert loc._3utr_indices == [(CDS_END, len(SEQ))]
    assert loc.start_codons == [(CDS_START, CDS_START + 3)]
    assert loc.stop_codons == [(CDS_END - 3, CDS_END)]
    assert len(CanonicalCds.from_locus(loc).sequence) == CDS_END - CDS_START


@pytest.mark.parametrize(
    "cds_start, cds_end",
    [
        (10, 39),  # not in-frame (29 nt)
        (40, 10),  # start >= end
        (-1, 30),  # out of range
        (10, 10),  # empty
        (10, None),  # only one bound given
        (None, 40),  # only one bound given
    ],
)
def test_invalid_cds_raises(cds_start, cds_end):
    with pytest.raises(ValueError):
        LocusInfo(seq=SEQ, cds_start=cds_start, cds_end=cds_end)
