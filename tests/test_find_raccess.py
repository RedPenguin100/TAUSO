import pytest

from tauso._raccess.core import find_raccess


def test_sanity():
    config = find_raccess()
    print(config)