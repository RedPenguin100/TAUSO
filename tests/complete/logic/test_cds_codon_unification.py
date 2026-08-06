"""The codon-usage CDS windows and the struct_sense_in_cds flag share one canonical
transcript model: a local CDS codon window is populated for exactly the ASOs flagged
struct_sense_in_cds=1.

The structure flags (canonical-transcript LocusInfo intervals) and the codon windows
must agree on which ASOs sit in the coding region; these tests pin that agreement.
"""

import pytest
import pandas as pd

from tauso.data.consts import STRUCT_SENSE_IN_CDS

LOCAL_COLS = [f"local_coding_region_around_ASO_{f}" for f in (20, 40, 70)]


# logic
def test_coding_region_pipeline_equivalence(final_data):
    """in_coding_region (codon pipeline) must equal struct_sense_in_cds (structure pipeline)."""
    struct_sense = final_data[STRUCT_SENSE_IN_CDS].fillna(0).astype(bool)
    in_coding = final_data["in_coding_region"].astype(bool)

    pd.testing.assert_series_equal(struct_sense, in_coding, check_names=False)


# logic
@pytest.mark.parametrize("local_coding_region", [f"local_coding_region_around_ASO_{f}" for f in (20, 40, 70)])
def test_window_presence_iff_in_coding(final_data, local_coding_region):
    """Window column is non-empty if and only if in_coding_region is True."""
    nonempty = final_data[local_coding_region].fillna("").str.len() > 0
    in_coding = final_data["in_coding_region"].astype(bool)

    pd.testing.assert_series_equal(nonempty, in_coding, check_names=False)

# logic
def test_cds_implies_coding_window(final_data):
    """Every struct_sense_in_cds=1 ASO must receive a local CDS window — guards the old
    longest-CDS path's regression that silently dropped real coding ASOs."""
    is_cds = final_data[STRUCT_SENSE_IN_CDS].fillna(False).astype(bool)
    windows = final_data["local_coding_region_around_ASO_40"].fillna("")

    # Isolate rows violating the rule: in CDS but window is empty
    mismatches = final_data[is_cds & (windows == "")]

    assert mismatches.empty, f"Missing windows for indices: {mismatches.index.tolist()}"

# logic
def test_local_window_length_in_frame(final_data):
    """Frame-aligned windows start on a codon boundary, so for the symmetric ±flank
    windows the populated CDS slice is consistent with the requested span."""
    windows = final_data["local_coding_region_around_ASO_40"].fillna("")
    populated = windows[windows.str.len() > 0]
    assert len(populated) > 0, "no populated CDS windows to check"
    # ±40 nt flank around a ~20 nt footprint, frame-aligned: comfortably under 120 nt.
    assert populated.str.len().max() <= 120
