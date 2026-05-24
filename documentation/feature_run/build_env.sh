#!/usr/bin/env bash
#
# Build the reproducible `tauso_repro` env from the explicit conda lock and
# editable-install tauso. Designed to be run on a SLURM compute node via:
#
#     tauso_run --cpu=8 --mem=32G bash <path-to-this-script>
#
# It is kept as a script (not an inline one-liner) because `tauso_run` does not
# pass `&&` chains / `|` pipes / `-n <env>` flags / quotes through cleanly.
#
# Override any of these via the environment if your layout differs:
#     MM                 path to the micromamba binary
#     MAMBA_ROOT_PREFIX  isolated root prefix for envs/packages
#     ENV                env name to create
set -euo pipefail

MM="${MM:-/tamir2/kovaliov/bin/micromamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/tamir2/kovaliov/micromamba}"
ENV="${ENV:-tauso_repro}"

# repo root = two levels up from this script (documentation/feature_run/..)
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO"

if ! "$MM" env list | grep -q "$ENV"; then
  echo ">>> creating env '$ENV' from conda-linux-64-dev.lock (explicit, pinned)"
  "$MM" create -y -n "$ENV" -f conda-linux-64-dev.lock
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
