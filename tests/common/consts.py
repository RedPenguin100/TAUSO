import os
from pathlib import Path

TESTS_PATH =  Path(os.path.dirname(os.path.dirname(__file__)))
UNIT_TESTS_PATH =  TESTS_PATH / 'unit'
UNIT_TESTS_CACHE = UNIT_TESTS_PATH / "cache"
INTEGRATION_TESTS_PATH = TESTS_PATH / 'integration'