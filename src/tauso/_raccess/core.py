import os
import shutil
from pathlib import Path

from ..data.data import get_data_dir
from .consts import RACCESS_EXE_ENV, RACCESS_EXE_NAME


def find_raccess():
    # 1. Environment variable override (direct path to executable)
    env = os.environ.get(RACCESS_EXE_ENV)
    if env and os.path.isfile(env) and os.access(env, os.X_OK):
        return env

    # 2. Look in the data directory (respects TAUSO_DATA_DIR)
    data_dir = Path(get_data_dir())
    bin_path = data_dir / "raccess" / "bin" / RACCESS_EXE_NAME
    if bin_path.exists() and os.access(bin_path, os.X_OK):
        return str(bin_path)

    # 3. PATH
    which = shutil.which(RACCESS_EXE_NAME)
    if which:
        return which

    raise RuntimeError(f"raccess not found. Run `tauso setup-raccess` or set {RACCESS_EXE_ENV} environment variable.")
