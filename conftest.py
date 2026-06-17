"""Root conftest: CLI options and marker gating.

These hooks must live at the rootdir. Since pytest 9.1, `pytest_addoption` is only
collected from the rootdir conftest (and the invocation/arg paths), not from conftests
in subdirectories. CI runs `pytest` from the repo root with no path argument, so an
`--n-jobs` option defined only in tests/conftest.py was never registered, failing with
"unrecognized arguments: --n-jobs".
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
