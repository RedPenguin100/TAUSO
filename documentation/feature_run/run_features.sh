#!/usr/bin/env bash
# Calculate features. Path-relative + self-locating. Usage:
#   run_features.sh            all features, single job
#   run_features.sh K N        gene partition K of N (saves to a partition dir)
#   run_features.sh merge N    merge N partition dirs into the main feature dir
# Env: MM, MAMBA_ROOT_PREFIX, TAUSO_DATA_DIR, ENV, CPUS (default nproc),
#      OFFHOME_CACHE=1 (keep caches off $HOME), OVERWRITE=1 (force recompute).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$(dirname "$REPO")"
export MAMBA_ROOT_PREFIX="$BASE/micromamba"
export TAUSO_DATA_DIR="$BASE/.tauso_data"
ENV="${ENV:-tauso_repro}"
if [ -z "${MM:-}" ]; then
  if [ -x "$BASE/bin/micromamba" ]; then MM="$BASE/bin/micromamba"; else MM="micromamba"; fi
fi
if [ "${OFFHOME_CACHE:-0}" = "1" ]; then
  export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$BASE/.cache}"
  export TMPDIR="${TMPDIR:-$BASE/.tmp}"
  mkdir -p "$XDG_CACHE_HOME" "$TMPDIR"
fi
CPUS="${CPUS:-$(nproc)}"
cd "$REPO"
# Fast guard: a -march=native raccess binary SIGILLs on a node unlike the build
# host. Catch it here (~1s) instead of after grinding to sense-accessibility.
# Skipped for `merge` (no raccess needed). Set SKIP_RACCESS_CHECK=1 to bypass.
RACCESS="$TAUSO_DATA_DIR/raccess/bin/run_raccess"
if [ "${1:-all}" != "merge" ] && [ "${SKIP_RACCESS_CHECK:-0}" != "1" ] && [ -x "$RACCESS" ]; then
  "$RACCESS" >/dev/null 2>&1 || if [ "$?" -eq 132 ]; then
    echo "ABORT: raccess SIGILL on $(hostname) -- rebuild portably before running:" >&2
    echo "  RACCESS_MARCH=x86-64-v2 tauso setup-raccess --force-clone --march x86-64-v2" >&2
    exit 1
  fi
fi
feats() { "$MM" run -n "$ENV" python -u -m notebooks.features.calculate_features --dataset oligo "$@"; }
case "${1:-all}" in
  merge) feats --merge "$2" ;;
  all)   feats --cpus "$CPUS" ${OVERWRITE:+--overwrite} ;;
  *)     feats --cpus "$CPUS" --partition "$1" "$2" ${OVERWRITE:+--overwrite} ;;
esac
