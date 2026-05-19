import os
import shutil
from pathlib import Path

from tauso._raccess.consts import RACCESS_EXE_ENV, RACCESS_EXE_NAME
from tauso.data.data import get_data_dir


def find_raccess():
    # 1. Environment variable override (direct path to executable)
    env = os.environ.get(RACCESS_EXE_ENV)
    if env and os.path.isfile(env) and os.access(env, os.X_OK):
        return env

    # 2. Look in the data directory (respects TAUSO_DATA_DIR)
    data_dir = Path(get_data_dir())
    data_path = data_dir / "raccess" / "src" / "raccess" / RACCESS_EXE_NAME
    if data_path.exists() and os.access(data_path, os.X_OK):
        return str(data_path)

    # 3. Look relative to this file (in case it's built in the repo/package)
    # e.g. tauso/_raccess/bin/run_raccess
    repo_path = Path(__file__).parent / "bin" / RACCESS_EXE_NAME
    if repo_path.exists() and os.access(repo_path, os.X_OK):
        return str(repo_path)

    # 4. PATH
    which = shutil.which(RACCESS_EXE_NAME)
    if which:
        return which

    raise RuntimeError(
        f"raccess not found. Run `tauso install-raccess` or set {RACCESS_EXE_ENV} environment variable."
    )
