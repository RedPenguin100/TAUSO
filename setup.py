import os
import subprocess
import shutil
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
from setuptools.command.egg_info import egg_info

# --- CONFIGURATION ---
PACKAGE_NAME = "tauso"
BINARY_NAME = "risearch_executable"

# ---------------------

def build_risearch_binary(force_target_dir=None):
    """
    Compiles the C binary and copies it.
    If force_target_dir is None, it defaults to 'src/tauso' (for editable).
    """
    root_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.join(root_dir, "external", "risearch", "RIsearch1")
    src_dir = os.path.join(base_dir, "src")
    bin_dir = os.path.join(base_dir, "bin")

    # 1. Determine Destination (Updated to include 'out')
    if force_target_dir:
        # Standard install (pip install .) -> goes to build/lib/tauso/out
        dest_dir = os.path.join(force_target_dir, "out")
    else:
        # Editable install (pip install -e .) -> goes to src/tauso/out
        dest_dir = os.path.join(root_dir, "src", PACKAGE_NAME, "out")

    dest_path = os.path.join(dest_dir, BINARY_NAME)

    # Optimization: If binary already exists in src during editable install, skip rebuild
    if not force_target_dir and os.path.exists(dest_path):
        print(f"Binary already exists at {dest_path}, skipping rebuild.")
        return

    print(f"--- BUILDING EXTENSION ---")

    # 2. Check Sources
    if not os.path.exists(src_dir):
        print(f"WARNING: C source not found at {src_dir}. Binary will not be built.")
        return

    # 3. Build
    os.makedirs(bin_dir, exist_ok=True)
    subprocess.check_call(["make", "-C", src_dir])

    # 4. Find Binary
    compiled_bin = os.path.join(bin_dir, "RIsearch")
    if not os.path.exists(compiled_bin):
        compiled_bin = os.path.join(src_dir, "RIsearch")

    if not os.path.exists(compiled_bin):
        raise FileNotFoundError("Make ran successfully, but RIsearch binary not found.")

    # 5. Copy
    # This creates src/tauso/out/ if it doesn't exist
    os.makedirs(dest_dir, exist_ok=True)
    print(f"Copying binary to {dest_path}")
    shutil.copy(compiled_bin, dest_path)
    os.chmod(dest_path, 0o755)


class CustomBuild(build_py):
    """Triggered by: pip install ."""

    def run(self):
        # Build into the temp build folder
        target_dir = os.path.join(self.build_lib, PACKAGE_NAME)
        build_risearch_binary(target_dir)
        super().run()


class CustomEggInfo(egg_info):
    """
    Triggered by BOTH: pip install . AND pip install -e .
    """

    def run(self):
        build_risearch_binary(force_target_dir=None)
        super().run()


setup(
    name="tauso",
    version="1.0",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    cmdclass={
        'build_py': CustomBuild,
        'egg_info': CustomEggInfo,
    },
    # Update package_data to include the 'out/' prefix
    package_data={PACKAGE_NAME: [f"out/{BINARY_NAME}"]},
    include_package_data=True,
)