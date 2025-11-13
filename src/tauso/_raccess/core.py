import os, shutil

from tauso._raccess.consts import RACCESS_EXE_ENV, RACCESS_EXE_NAME, TAUSO_SHARE


def find_raccess():
    env = os.environ.get(RACCESS_EXE_ENV)
    if env and os.path.isfile(env) and os.access(env, os.X_OK):
        return env

    default = TAUSO_SHARE / 'raccess' / 'src' /'raccess' / RACCESS_EXE_NAME
    if default.exists() and os.access(default, os.X_OK):
        return str(default)

    which = shutil.which(RACCESS_EXE_NAME)
    if which:
        return which

    raise RuntimeError("raccess not found. Run `tauso-setup-raccess` or set RACCESS_EXE.")