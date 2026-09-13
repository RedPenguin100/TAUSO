"""The model cache key has to follow the registry, or a stale cache answers for the wrong model."""

import importlib

from click.testing import CliRunner

from tauso import cli
from tauso.inference.predict import MODEL_FILES, ZENODO_MODEL_RECORD, model_cache_key


def test_the_key_names_the_record_it_fetches_from():
    assert ZENODO_MODEL_RECORD in model_cache_key()


def test_the_key_moves_when_a_model_does(monkeypatch):
    before = model_cache_key()
    moved = {name: dict(spec) for name, spec in MODEL_FILES.items()}
    moved["v1"] = {"filename": "tauso_score_v9.json", "md5": "0" * 32}
    # By name it resolves to the function tauso.inference re-exports, not the module.
    monkeypatch.setattr(importlib.import_module("tauso.inference.predict"), "MODEL_FILES", moved)

    assert model_cache_key() != before


def test_the_key_holds_still_when_nothing_moves():
    assert model_cache_key() == model_cache_key()


def test_the_command_prints_it():
    result = CliRunner().invoke(cli.main, ["model-cache-key"])

    assert result.exit_code == 0
    assert result.output.strip() == model_cache_key()
