"""Root conftest: CLI options and marker gating.

pytest collects `pytest_addoption` only from the rootdir conftest and the paths named on
the command line, not from subdirectory conftests, so these option hooks must live here at
the repo root -- the test command runs `pytest` from the root with no path argument.
"""

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run-integration",
        action="store_true",
        default=False,
        help="Run integration tests",
    )
    parser.addoption(
        "--run-benchmarks",
        action="store_true",
        default=False,
        help="Run benchmark tests",
    )
    parser.addoption(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of parallel workers passed to populate functions (default: 1)",
    )


def pytest_runtest_setup(item):
    if "integration" in item.keywords and not item.config.getoption("--run-integration"):
        pytest.skip("Skipping integration tests (use --run-integration to run)")
    if "benchmark" in item.keywords and not item.config.getoption("--run-benchmarks"):
        pytest.skip("Skipping benchmark tests (use --run-benchmarks to run)")
