"""Per-ASO fold template selection: pre-mRNA (default) vs canonical CDS, chosen per row."""

import pandas as pd
import pytest

from tauso.data.consts import CANONICAL_GENE_NAME, STRUCTURE_SENSE_LENGTH, STRUCTURE_SENSE_START
from tauso.genome.LocusInfo import LocusInfo
from tauso.populate.populate_fold import (
    SEQUENCE_SOURCE_CANONICAL_CDS,
    SEQUENCE_SOURCE_PREMRNA,
    _build_templates_by_source,
    _lightweight_gene_to_data,
    _select_source,
    populate_mfe_features,
)

SEQ = "ACGTACGT" * 20  # 160 nt, single-exon custom transcript
CDS_START, CDS_END = 10, 40  # 30 nt in-frame coding span
CDS_SEQ = SEQ[CDS_START:CDS_END]

SOURCE_COL = "fold_source"


def _coding_locus(seq=SEQ, cds_start=CDS_START, cds_end=CDS_END):
    return LocusInfo(seq=seq, cds_start=cds_start, cds_end=cds_end)


# --- template builder (per source) -------------------------------------------------


def test_premrna_returns_full_mrna():
    out = _lightweight_gene_to_data(["G"], {"G": _coding_locus()}, SEQUENCE_SOURCE_PREMRNA)
    assert out == {"G": SEQ}


def test_canonical_cds_returns_spliced_cds():
    out = _lightweight_gene_to_data(["G"], {"G": _coding_locus()}, SEQUENCE_SOURCE_CANONICAL_CDS)
    assert out == {"G": CDS_SEQ}


def test_canonical_cds_omits_non_coding():
    out = _lightweight_gene_to_data(["G"], {"G": LocusInfo(seq=SEQ)}, SEQUENCE_SOURCE_CANONICAL_CDS)
    assert out == {}


def test_unknown_source_raises():
    with pytest.raises(ValueError):
        _lightweight_gene_to_data(["G"], {"G": _coding_locus()}, "spliced_mrna")


# --- per-row selection logic -------------------------------------------------------


def test_build_templates_default_is_premrna_only():
    df = pd.DataFrame({SOURCE_COL: [SEQUENCE_SOURCE_CANONICAL_CDS]})
    tbs = _build_templates_by_source(["G"], {"G": _coding_locus()}, df, sequence_source_col=None)
    assert set(tbs) == {SEQUENCE_SOURCE_PREMRNA}


def test_build_templates_covers_column_values():
    df = pd.DataFrame({SOURCE_COL: [SEQUENCE_SOURCE_CANONICAL_CDS, None]})
    tbs = _build_templates_by_source(["G"], {"G": _coding_locus()}, df, sequence_source_col=SOURCE_COL)
    assert set(tbs) == {SEQUENCE_SOURCE_PREMRNA, SEQUENCE_SOURCE_CANONICAL_CDS}
    assert tbs[SEQUENCE_SOURCE_PREMRNA] == {"G": SEQ}
    assert tbs[SEQUENCE_SOURCE_CANONICAL_CDS] == {"G": CDS_SEQ}


@pytest.mark.parametrize(
    "value, expected",
    [
        (SEQUENCE_SOURCE_CANONICAL_CDS, SEQUENCE_SOURCE_CANONICAL_CDS),
        (SEQUENCE_SOURCE_PREMRNA, SEQUENCE_SOURCE_PREMRNA),
        (None, SEQUENCE_SOURCE_PREMRNA),  # blank -> fallback
        ("nonsense", SEQUENCE_SOURCE_PREMRNA),  # unrecognized -> fallback
    ],
)
def test_select_source(value, expected):
    tbs = {SEQUENCE_SOURCE_PREMRNA: {}, SEQUENCE_SOURCE_CANONICAL_CDS: {}}
    assert _select_source({SOURCE_COL: value}, tbs, SOURCE_COL) == expected


def test_select_source_none_column_is_premrna():
    assert _select_source({}, {SEQUENCE_SOURCE_PREMRNA: {}}, None) == SEQUENCE_SOURCE_PREMRNA


# --- end-to-end: the template is chosen per ASO within a single call ----------------

MFE_COL = "fold_mfe_win25_flank30_step4"
SETTINGS = [(30, 25, 4)]
# A start valid on the 160-nt pre-mRNA but past the end of the 30-nt canonical CDS, so
# routing shows up as folded-vs-NaN rather than a structure-dependent MFE difference.
OUT_OF_CDS_START = 100


def _row(source):
    return {
        CANONICAL_GENE_NAME: "G",
        STRUCTURE_SENSE_START: OUT_OF_CDS_START,
        STRUCTURE_SENSE_LENGTH: 20,
        SOURCE_COL: source,
    }


def test_mfe_template_chosen_per_row():
    gene_to_data = {"G": _coding_locus()}  # pre-mRNA 160 nt, canonical CDS 30 nt
    df = pd.DataFrame([_row(SEQUENCE_SOURCE_PREMRNA), _row(SEQUENCE_SOURCE_CANONICAL_CDS)])
    out, _ = populate_mfe_features(df, gene_to_data, settings=SETTINGS, sequence_source_col=SOURCE_COL)
    # Identical inputs; only the per-row source differs -> only the CDS-routed row runs out of template.
    assert pd.notna(out.loc[0, MFE_COL])
    assert pd.isna(out.loc[1, MFE_COL])


def test_mfe_default_ignores_source_column():
    gene_to_data = {"G": _coding_locus()}
    df = pd.DataFrame([_row(SEQUENCE_SOURCE_PREMRNA), _row(SEQUENCE_SOURCE_CANONICAL_CDS)])
    out, _ = populate_mfe_features(df, gene_to_data, settings=SETTINGS)  # no sequence_source_col
    # Without routing, both rows fold the pre-mRNA at the same start -> both non-NaN and identical.
    assert pd.notna(out.loc[0, MFE_COL])
    assert out.loc[0, MFE_COL] == out.loc[1, MFE_COL]
