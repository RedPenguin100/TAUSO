import pytest

from tauso.target_finder import get_3utr_gfp

"""
File for tests to check that we got the right sequences
"""


def test_gfp_3_utr():
    gfp_3_utr = get_3utr_gfp()

    assert gfp_3_utr is not None
    assert gfp_3_utr[0:3] == "TGG"
    assert gfp_3_utr[-3:] == "GCA"
