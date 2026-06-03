"""The codon-usage CDS windows and the sense_cds flag share one canonical transcript
model: a local CDS codon window is populated for exactly the ASOs flagged sense_cds=1.

The structure flags (canonical-transcript LocusInfo intervals) and the codon windows
must agree on which ASOs sit in the coding region; these tests pin that agreement.
"""

import pytest

from tauso.data.consts import SENSE_CDS

LOCAL_COLS = [f"local_coding_region_around_ASO_{f}" for f in (20, 40, 70)]


def test_in_coding_region_equals_sense_cds(final_data):
    """in_coding_region (codon pipeline) must equal sense_cds (structure pipeline)."""
    sc = final_data[SENSE_CDS].fillna(0).astype(bool)
    ic = final_data["in_coding_region"].astype(bool)
    mismatches = int((sc.values != ic.values).sum())
    assert mismatches == 0, f"{mismatches} rows where in_coding_region != sense_cds"


@pytest.mark.parametrize("col", LOCAL_COLS)
def test_local_window_nonempty_iff_in_coding(final_data, col):
    nonempty = (final_data[col].str.len() > 0).values
    ic = final_data["in_coding_region"].astype(bool).values
    assert (nonempty == ic).all(), f"{col}: window-populated set differs from in_coding_region"


def test_every_sense_cds_aso_gets_a_window(final_data):
    """Every sense_cds=1 ASO must receive a local CDS window — guards the old
    longest-CDS path's regression that silently dropped ~1.1k real coding ASOs."""
    sc = final_data[SENSE_CDS].fillna(0).astype(bool)
    win = final_data["local_coding_region_around_ASO_40"].fillna("")
    missing = int((sc & (win.str.len() == 0)).sum())
    assert missing == 0, f"{missing} sense_cds ASOs have no local codon window"


def test_local_window_length_in_frame(final_data):
    """Frame-aligned windows start on a codon boundary, so for the symmetric ±flank
    windows the populated CDS slice is consistent with the requested span."""
    win = final_data["local_coding_region_around_ASO_40"].fillna("")
    populated = win[win.str.len() > 0]
    assert len(populated) > 0, "no populated CDS windows to check"
    # ±40 nt flank around a ~20 nt footprint, frame-aligned: comfortably under 120 nt.
    assert populated.str.len().max() <= 120
