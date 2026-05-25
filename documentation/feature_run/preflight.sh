#!/usr/bin/env bash
#
# Preflight for the feature run: fail fast on missing prerequisites BEFORE
# launching expensive multi-node jobs. Path-relative + self-locating.
#
# Run as ONE cheap job, then only launch the 4x128 partitions if it passes:
#     tauso_run --cpu=16 "bash <BASE>/TAUSO/documentation/feature_run/preflight.sh"
#
# It checks two things that have bitten us:
#   1. raccess runs on THIS node (a -march=native binary SIGILLs on a node
#      unlike the build host) -- caught in ~1s.
#   2. a 1-gene run through the FULL pipeline succeeds -- exercises every
#      calculator and every data file (mrna_half_life.parquet, genome, depmap,
#      ...) via the real code path, so any missing prerequisite shows up in
#      ~1-2 min instead of after a multi-node launch.
#
# Env overrides: MM, MAMBA_ROOT_PREFIX, TAUSO_DATA_DIR, ENV, CPUS, CANARY_N.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE="$(dirname "$REPO")"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$BASE/micromamba}"
export TAUSO_DATA_DIR="${TAUSO_DATA_DIR:-$BASE/.tauso_data}"
ENV="${ENV:-tauso_repro}"
if [ -z "${MM:-}" ]; then
  if [ -x "$BASE/bin/micromamba" ]; then MM="$BASE/bin/micromamba"; else MM="micromamba"; fi
fi
CPUS="${CPUS:-$(nproc)}"
cd "$REPO"
run() { "$MM" run -n "$ENV" "$@"; }

echo ">>> preflight: TAUSO_DATA_DIR=$TAUSO_DATA_DIR  (node $(hostname))"

# --- 1. raccess: exists and actually runs on this node (catches SIGILL) -------
RACCESS="$TAUSO_DATA_DIR/raccess/bin/run_raccess"
if [ ! -x "$RACCESS" ]; then
  echo "FAIL: raccess binary missing at $RACCESS"
  echo "      Build it:  RACCESS_MARCH=x86-64-v2 tauso setup-raccess --force-clone --march x86-64-v2"
  exit 1
fi
"$RACCESS" >/dev/null 2>&1; rc=$?
if [ "$rc" -eq 132 ]; then   # 128 + SIGILL(4)
  echo "FAIL: raccess SIGILL on $(hostname) -- binary built for a different CPU."
  echo "      Rebuild portably (note --force-clone, else the build short-circuits):"
  echo "      RACCESS_MARCH=x86-64-v2 tauso setup-raccess --force-clone --march x86-64-v2"
  exit 1
fi
echo "  ok    raccess runs on this node (exit $rc, no SIGILL)"

# --- 2. 1-gene canary through the full pipeline -------------------------------
# genes[0::CANARY_N] is exactly one gene when CANARY_N exceeds the gene count.
CANARY_N="${CANARY_N:-100000}"
CANARY_DIR="$REPO/notebooks/saved_features_oligo_p0of${CANARY_N}"
rm -rf "$CANARY_DIR"
echo ">>> canary: 1-gene full pipeline (partition 0 of $CANARY_N), cpus=$CPUS ..."
if run python -u -m notebooks.features.calculate_features \
      --dataset oligo --cpus "$CPUS" --partition 0 "$CANARY_N" --overwrite; then
  rm -rf "$CANARY_DIR"
  echo ">>> PREFLIGHT PASSED -- safe to launch the partition jobs."
else
  rc=$?
  rm -rf "$CANARY_DIR"
  echo ">>> PREFLIGHT FAILED (exit $rc) -- fix the error above before launching."
  exit 1
fi
