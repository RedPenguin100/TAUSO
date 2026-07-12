"""Terminal one-hot encoding: 5' and length-aware 3' windows."""

import pandas as pd

from tauso.data.consts import ASO_SEQUENCE
from tauso.populate.populate_sequence import (
    TERMINAL_OHE_N,
    one_hot_feature_names,
    populate_sequence_one_hot_encoded,
)

N = TERMINAL_OHE_N


def _encode(seqs):
    df = pd.DataFrame({ASO_SEQUENCE: seqs})
    out, names = populate_sequence_one_hot_encoded(df)
    return out, names


def test_feature_names_shape():
    names = one_hot_feature_names()
    assert len(names) == 2 * N * 4
    # 5' block first, then 3' block; four nucleotide columns per position
    assert names[:4] == [f"ohe_pos0_{n}" for n in "ACGT"]
    assert names[N * 4:N * 4 + 4] == [f"ohe_3p0_{n}" for n in "ACGT"]
    assert set(names) == {f"ohe_pos{p}_{n}" for p in range(N) for n in "ACGT"} | {
        f"ohe_3p{p}_{n}" for p in range(N) for n in "ACGT"
    }


def test_populate_returns_expected_columns():
    out, names = _encode(["ACGTACGTACGTACGT"])
    assert names == one_hot_feature_names()
    assert all(c in out.columns for c in names)


def test_five_prime_window():
    # 5' positions read left-to-right from the start
    out, _ = _encode(["ACGTACGTAC"])  # 10 nt, longer than N
    assert out.loc[0, "ohe_pos0_A"] == 1 and out.loc[0, "ohe_pos0_C"] == 0
    assert out.loc[0, "ohe_pos1_C"] == 1
    assert out.loc[0, "ohe_pos2_G"] == 1
    assert out.loc[0, "ohe_pos3_T"] == 1


def test_three_prime_window_is_length_aware():
    # Two sequences with the SAME 3' end but different lengths must share ohe_3p*.
    a = "AAAAAAAAAA" + "TGCA"  # 3' end ...T G C A
    b = "GG" + "TGCA"          # same 3' end, much shorter
    out, _ = _encode([a, b])
    for col in ("ohe_3p0_A", "ohe_3p1_C", "ohe_3p2_G", "ohe_3p3_T"):
        assert out.loc[0, col] == 1, col
        assert out.loc[0, col] == out.loc[1, col], col
    # ohe_3p0 is the very last nucleotide (A here), not a 5'-indexed position
    assert out.loc[0, "ohe_3p0_A"] == 1
    assert out.loc[0, "ohe_3p1_C"] == 1  # second-from-last


def test_short_sequence_is_zero_padded_not_erroring():
    # Sequence shorter than N: out-of-range terminal slots are all-zero, no crash.
    out, _ = _encode(["ACG"])  # 3 nt < N=6
    assert out.loc[0, "ohe_pos0_A"] == 1
    # position 5 from the 5' end does not exist -> all four channels zero
    assert out.loc[0, [f"ohe_pos5_{n}" for n in "ACGT"]].sum() == 0
    # 3' end: ohe_3p0 is G (last nt); slots beyond the sequence are zero
    assert out.loc[0, "ohe_3p0_G"] == 1
    assert out.loc[0, [f"ohe_3p5_{n}" for n in "ACGT"]].sum() == 0


def test_u_maps_like_t():
    out, _ = _encode(["ACGUACGU"])
    # U at 5' position 3 encodes on the T channel
    assert out.loc[0, "ohe_pos3_T"] == 1
