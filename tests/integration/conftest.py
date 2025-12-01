import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--run-integration",
        action="store_true",
        default=False,
        help="Run integration tests"
    )

def pytest_runtest_setup(item):
    if "integration" in item.keywords and not item.config.getoption("--run-integration"):
        pytest.skip("Skipping integration tests (use --run-integration to run)")
