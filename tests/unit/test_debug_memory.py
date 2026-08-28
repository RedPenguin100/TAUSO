import logging

import numpy as np
import pandas as pd
import pytest

from tauso.data.consts import CANONICAL_GENE_NAME, CHEMICAL_PATTERN, PS_PATTERN
from tauso.debug import log_dataframe_memory
from tauso.populate.calculators.calculator import Calculator


@pytest.fixture
def frame():
    return pd.DataFrame(
        {
            "index_oligo": np.arange(100),
            "aso_sequence": ["ACGT" * 5] * 100,
            CANONICAL_GENE_NAME: ["MALAT1"] * 100,
            PS_PATTERN: ["sssss"] * 100,
            CHEMICAL_PATTERN: ["eeeddd"] * 100,
        }
    )


def test_log_dataframe_memory_reports_total_and_heaviest_columns(frame, caplog):
    with caplog.at_level(logging.INFO, logger="tauso.debug"):
        log_dataframe_memory(frame, "after calculate_demo")

    message = caplog.text
    assert "after calculate_demo" in message
    assert f"over {frame.shape[1]} columns" in message
    # the string column outweighs the integer one, so it must lead the heaviest list
    assert message.index("aso_sequence") < message.index("index_oligo")


def _stub_pipeline(calc, monkeypatch):
    """Replace every calculate_* step with a no-op so calculate_all is cheap."""
    for name in dir(calc):
        if name.startswith("calculate_") and name != "calculate_all":
            monkeypatch.setattr(calc, name, lambda: None)


def test_calculate_all_profiles_memory_only_when_enabled(frame, monkeypatch, caplog):
    monkeypatch.delenv("TAUSO_PROFILE_MEM", raising=False)
    calc = Calculator(data=frame, data_version=None, overwrite=True)
    _stub_pipeline(calc, monkeypatch)

    with caplog.at_level(logging.INFO, logger="tauso.debug"):
        calc.calculate_all()
    assert "[MEM]" not in caplog.text

    monkeypatch.setenv("TAUSO_PROFILE_MEM", "1")
    caplog.clear()
    with caplog.at_level(logging.INFO, logger="tauso.debug"):
        calc.calculate_all()
    assert caplog.text.count("[MEM]") > 1
