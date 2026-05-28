#!/usr/bin/env bash
#
# Build the reproducible tauso env from the explicit conda lock and
# editable-install tauso.
#
# Path-relative and self-locating: the repo is the clone two levels above this
# script (<BASE>/TAUSO/documentation/feature_run/build_env.sh), so <BASE> is the
# directory above the clone. Nothing is hard-coded; launch it however you like
# (directly, or wrapped by a job scheduler):
#
#     bash <BASE>/TAUSO/documentation/feature_run/build_env.sh
#
# Override via the environment if your layout differs:
#     MM                 micromamba binary (default: $BASE/bin/micromamba, else PATH)
#     MAMBA_ROOT_PREFIX  env/package store (default: $BASE/micromamba)
#     ENV                env name           (default: tauso_repro)
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # the TAUSO clone
BASE="$(dirname "$REPO")"                                    # the chosen base dir

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$BASE/micromamba}"
ENV="${ENV:-tauso_repro}"
if [ -z "${MM:-}" ]; then
  if [ -x "$BASE/bin/micromamba" ]; then MM="$BASE/bin/micromamba"; else MM="micromamba"; fi
fi

# Optional: on machines with a quota-limited $HOME (e.g. some clusters), uv/pip
# default to ~/.cache and can blow the home quota mid-install. Set OFFHOME_CACHE=1
# to keep all caches/temp under the base dir instead. Off by default.
if [ "${OFFHOME_CACHE:-0}" = "1" ]; then
  export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$BASE/.cache}"
  export UV_CACHE_DIR="${UV_CACHE_DIR:-$BASE/.uv_cache}"
  export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$BASE/.pip_cache}"
  export TMPDIR="${TMPDIR:-$BASE/.tmp}"
  mkdir -p "$XDG_CACHE_HOME" "$UV_CACHE_DIR" "$PIP_CACHE_DIR" "$TMPDIR"
fi

cd "$REPO"

# Check by env prefix dir, NOT `env list | grep $ENV` — that one is fooled when $ENV
# is a substring of $MAMBA_ROOT_PREFIX's path (env list prints the root prefix path),
# which silently skips creation and later fails with "given prefix does not exist".
if [ ! -d "$MAMBA_ROOT_PREFIX/envs/$ENV" ]; then
  echo ">>> creating env '$ENV' from conda-linux-64-dev.lock (explicit, pinned)"
  "$MM" create -y -n "$ENV" --file conda-linux-64-dev.lock
else
  echo ">>> env '$ENV' already exists; skipping create"
fi

echo ">>> editable install of tauso into '$ENV'"
if "$MM" run -n "$ENV" uv --version >/dev/null 2>&1; then
  "$MM" run -n "$ENV" uv pip install -e ".[dev]"
else
  "$MM" run -n "$ENV" pip install -e ".[dev]"
fi

echo ">>> verify"
"$MM" run -n "$ENV" tauso --help | head
"$MM" run -n "$ENV" python -c "import tauso; print('tauso OK', tauso.__file__)"
