import os
import subprocess
import pytest
from tauso._raccess.core import find_raccess
from tauso.data.data import get_data_dir


@pytest.fixture
def temp_data_dir(tmp_path):
    """Fixture to provide a clean TAUSO_DATA_DIR."""
    original_env = os.environ.get("TAUSO_DATA_DIR")
    os.environ["TAUSO_DATA_DIR"] = str(tmp_path)
    yield tmp_path
    if original_env:
        os.environ["TAUSO_DATA_DIR"] = original_env
    else:
        del os.environ["TAUSO_DATA_DIR"]


@pytest.mark.integration
def test_setup_raccess_custom_dir(temp_data_dir):
    """Test that setup-raccess builds raccess into TAUSO_DATA_DIR."""
    # Run the setup command
    # We use subprocess to simulate a real CLI call
    cmd = ["tauso", "setup-raccess"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)

    # 1. Check it exists in the data dir
    raccess_exe = temp_data_dir / "raccess" / "bin" / "run_raccess"
    assert raccess_exe.exists()
    assert os.access(raccess_exe, os.X_OK)

    # 2. Verify find_raccess identifies it. raccess is never bundled in the
    # package (non-distribution license); it lives only in TAUSO_DATA_DIR.
    found_path = find_raccess()
    assert found_path == str(raccess_exe)


def test_find_raccess_env_priority(temp_data_dir, tmp_path):
    """Test that RACCESS_EXE env var has the highest priority."""
    # Create a dummy executable elsewhere
    ext_dir = tmp_path / "external"
    ext_dir.mkdir()
    ext_exe = ext_dir / "my_raccess"
    ext_exe.write_text("#!/bin/sh\necho 'external'")
    ext_exe.chmod(0o755)

    os.environ["RACCESS_EXE"] = str(ext_exe)
    try:
        found_path = find_raccess()
        assert found_path == str(ext_exe)
    finally:
        del os.environ["RACCESS_EXE"]
