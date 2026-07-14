"""Locating splice-junction ASOs on the canonical CDS and recomputing their fold there."""

import numpy as np
import pandas as pd
import pytest

from tauso.data.consts import CANONICAL_GENE_NAME, STRUCTURE_SENSE_START
from tauso.genome.LocusInfo import LocusInfo
from tauso.populate.junction_fold import FOLD_SOURCE_COL, locate_on_cds, recompute_fold_on_cds
from tauso.util import get_antisense_rna


def _aperiodic(n):
    bases = "ACGT"
    return "".join(bases[(i * 5 + i // 3 + (i % 7)) % 4] for i in range(n))


SEQ = _aperiodic(160)
CDS_START, CDS_END = 12, 42  # 30 nt in-frame coding span
CDS = SEQ[CDS_START:CDS_END]
GENE = "G"


def _gene_to_data():
    return {GENE: LocusInfo(seq=SEQ, cds_start=CDS_START, cds_end=CDS_END)}


def _aso_targeting_cds(pos, length):
    """An ASO whose antisense matches CDS[pos:pos+length] (so it locates at CDS-frame `pos`)."""
    target_rna = CDS[pos : pos + length].upper().replace("T", "U")
    return get_antisense_rna(target_rna)


def test_locate_on_cds_returns_cds_frame_start():
    aso = _aso_targeting_cds(pos=5, length=20)
    starts = locate_on_cds([aso], [GENE], _gene_to_data())
    assert list(starts) == [5]


def test_locate_on_cds_absent_target_is_minus_one():
    starts = locate_on_cds(["A" * 20], [GENE], _gene_to_data())  # not present in CDS
    assert list(starts) == [-1]


def test_locate_on_cds_unknown_gene_and_noncoding():
    gtd = {GENE: LocusInfo(seq=SEQ)}  # non-coding -> empty CDS
    aso = _aso_targeting_cds(pos=5, length=20)
    assert list(locate_on_cds([aso], [GENE], gtd)) == [-1]  # empty CDS
    assert list(locate_on_cds([aso], ["ABSENT"], _gene_to_data())) == [-1]  # gene missing


def test_locate_groups_multiple_asos_per_gene():
    asos = [_aso_targeting_cds(3, 18), _aso_targeting_cds(9, 18), "T" * 18]
    starts = locate_on_cds(asos, [GENE, GENE, GENE], _gene_to_data())
    assert list(starts) == [3, 9, -1]


def test_recompute_fold_on_cds_fills_located_rows():
    df = pd.DataFrame(
        {
            "aso_sequence": [_aso_targeting_cds(4, 20), "A" * 20],  # located, not-located
            CANONICAL_GENE_NAME: [GENE, GENE],
        }
    )
    out, names = recompute_fold_on_cds(
        df, _gene_to_data(), sequence_col="aso_sequence", with_access=False, mfe_settings=[(30, 25, 4)]
    )
    assert names == ["fold_mfe_win25_flank30_step4"]
    assert out.loc[0, STRUCTURE_SENSE_START] == 4  # CDS-frame start
    assert out.loc[1, STRUCTURE_SENSE_START] == -1
    assert (out[FOLD_SOURCE_COL] == "canonical_cds").all()
    assert pd.notna(out.loc[0, names[0]])  # located row folds on the CDS
    assert pd.isna(out.loc[1, names[0]])  # unlocated row stays NaN
