"""Shape checks for every Calculator step.

`tests/complete` exercises the populate_* functions directly, so the Calculator's own
wiring — the expected-feature lists, the lazy imports, the skip path — has almost no
coverage. These tests drive every step with nothing missing, which runs each method's
expected-feature construction and its imports without doing any real computation.
"""

import logging

import numpy as np
import pandas as pd
import pytest

from tauso.data.consts import CANONICAL_GENE_NAME, CHEMICAL_PATTERN, PS_PATTERN
from tauso.populate.calculators.calculator import Calculator

# calculate_oligowalk has no missing-features check by design (it always recomputes) and
# shells out to the OligoWalk binary, so it cannot be driven this way.
STEPS = sorted(
    name
    for name in dir(Calculator)
    if name.startswith("calculate_") and name not in ("calculate_all", "calculate_oligowalk")
)


@pytest.fixture
def calc():
    data = pd.DataFrame(
        {
            "index_oligo": np.arange(10),
            CANONICAL_GENE_NAME: ["MALAT1"] * 10,
            PS_PATTERN: ["sssss"] * 10,
            CHEMICAL_PATTERN: ["eeeddd"] * 10,
        }
    )
    return Calculator(data=data, data_version=None, overwrite=False)


def test_every_step_is_discovered():
    # Guards against the list silently emptying if the naming convention changes.
    assert len(STEPS) >= 25


@pytest.mark.parametrize("step_name", STEPS)
def test_step_skips_cleanly_when_nothing_is_missing(calc, step_name, monkeypatch, caplog):
    """With every expected feature already on disk, a step must skip without computing.

    This still evaluates each step's expected-feature list and its lazy imports, which is
    where a refactor is most likely to break something.
    """
    monkeypatch.setattr(calc, "_get_missing_features", lambda expected: [])
    # Nothing should reach the populate layer; a step that tries is doing work it shouldn't.
    monkeypatch.setattr(calc, "_save_calculated_feature", _fail_on_save)
    # calculate_structure reloads its columns on the skip path; that needs a real feature
    # store, which is not what this test is about.
    monkeypatch.setattr(calc, "_load_features_into_data", lambda names: None)

    with caplog.at_level(logging.INFO):
        getattr(calc, step_name)()

    # calculate_structure says "Loading from disk" instead of "Skipping"; both report
    # that the features already exist, which is what matters here.
    assert "exist." in caplog.text


def _fail_on_save(*args, **kwargs):
    raise AssertionError("a step with nothing missing must not save features")
