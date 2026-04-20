import pytest
import subprocess

from notebooks.competitors.pfred_utils import populate_pfred, validate_docker_container


@pytest.fixture(scope="session")
def pfred_container():
    """
    Ensures the PFRED docker container is running before executing tests.
    If it isn't running, tests requiring it will be skipped gracefully.
    """
    container_name = "pfred"

    # Run the validation silently for tests
    is_running = validate_docker_container(container_name=container_name, verbose=False)

    if is_running:
        return container_name
    else:
        # pytest.skip aborts the test gracefully rather than failing it
        pytest.skip(
            f"Docker container '{container_name}' is not running, stopped, or missing. Skipping PFRED integration tests."
        )


@pytest.mark.skip(reason="TODO: fix this in the CI, currently too much setup with the docker probably")
@pytest.mark.parametrize("mini_sampled_data", [10000], indirect=True)
def test_pfred_features_regression(mini_sampled_data, pfred_container, dataframe_regression, tmp_path):
    """
    Tests the PFRED batch population against known regression baselines.
    Requires the pfred_container fixture to ensure Docker is running.
    """
    data = mini_sampled_data.copy()

    # Run the population step, passing tmp_path so file IO is isolated
    processed_data, feature_cols = populate_pfred(data=data, container_name=pfred_container, temp_dir=str(tmp_path))

    # Check only the newly generated feature columns against the baseline
    dataframe_regression.check(processed_data[feature_cols])
