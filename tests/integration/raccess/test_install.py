import os
import shutil
import subprocess
from pathlib import Path
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


@pytest.fixture
def clean_bin_dir():
    """Fixture to ensure the package bin directory is clean before/after tests."""
    bin_dir = Path(__file__).parents[3] / "src" / "tauso" / "_raccess" / "bin"
    temp_bin = None
    if bin_dir.exists():
        temp_bin = bin_dir.with_suffix(".bak")
        shutil.move(bin_dir, temp_bin)

    yield bin_dir

    if bin_dir.exists():
        shutil.rmtree(bin_dir)
    if temp_bin and temp_bin.exists():
        shutil.move(temp_bin, bin_dir)


@pytest.mark.integration
def test_setup_raccess_custom_dir(temp_data_dir, clean_bin_dir):
    """Test that setup-raccess installs into TAUSO_DATA_DIR and copies to bin/."""
    # Run the setup command
    # We use subprocess to simulate a real CLI call
    cmd = ["tauso", "setup-raccess"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)

    # 1. Check it exists in the data dir
    raccess_exe = temp_data_dir / "raccess" / "bin" / "run_raccess"
    assert raccess_exe.exists()
    assert os.access(raccess_exe, os.X_OK)

    # 2. Check it was copied to the package bin dir
    package_bin_exe = clean_bin_dir / "run_raccess"
    assert package_bin_exe.exists()
    assert os.access(package_bin_exe, os.X_OK)

    # 3. Verify find_raccess identifies it
    found_path = find_raccess()
    # It should find the one in the data dir first based on my implementation
    assert found_path == str(raccess_exe)


def test_find_raccess_fallback_to_bin(temp_data_dir, clean_bin_dir):
    """Test that find_raccess falls back to the package bin if data dir is empty."""
    # Create a dummy executable in the bin dir
    clean_bin_dir.mkdir(parents=True, exist_ok=True)
    dummy_exe = clean_bin_dir / "run_raccess"
    dummy_exe.write_text("#!/bin/sh\necho 'dummy'")
    dummy_exe.chmod(0o755)

    # Ensure data dir doesn't have it
    raccess_in_data = temp_data_dir / "raccess" / "bin" / "run_raccess"
    assert not raccess_in_data.exists()

    # Should find the one in bin/
    found_path = find_raccess()
    assert found_path == str(dummy_exe)


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
