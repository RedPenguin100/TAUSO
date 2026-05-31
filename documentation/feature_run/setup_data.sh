#!/usr/bin/env bash
#
# Data setup: build every reference dataset the feature pipeline needs, in the
# correct, non-redundant order. Path-relative and self-locating; launch however
# you like (directly, or wrapped by a job scheduler).
#
# Verified against tauso cli.py: `setup-all` already chains
#   genome -> bowtie -> omics(depmap, mrna-halflife, tgcn, ATtRACT, ribo-seq)
#         -> raccess
# so the omics sub-steps must NOT be invoked separately. The only repeat is
# `setup-tgcn` (also run by build-cell-context); it is idempotent.
#
# Overrides (env):
#   MM, MAMBA_ROOT_PREFIX, ENV, TAUSO_DATA_DIR
#   BOWTIE_THREADS         (default: nproc — the cores available to the process)
#   BOWTIE_MEM_PER_THREAD  (default: 2048 MB)
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # the TAUSO clone
BASE="$(dirname "$REPO")"                                    # the chosen base dir

export MAMBA_ROOT_PREFIX="$BASE/micromamba"
export TAUSO_DATA_DIR="$BASE/.tauso_data"
ENV="${ENV:-tauso_repro}"
if [ -z "${MM:-}" ]; then
  if [ -x "$BASE/bin/micromamba" ]; then MM="$BASE/bin/micromamba"; else MM="micromamba"; fi
fi
THREADS="${BOWTIE_THREADS:-$(nproc)}"
MEM_PER_THREAD="${BOWTIE_MEM_PER_THREAD:-2048}"

# raccess is compiled by setup-all. -march=native produces a binary that only
# runs on a CPU like the build node's, which SIGILLs on a heterogeneous cluster
# where build and run nodes differ. Default to the portable x86-64-v2 baseline
# (runs on essentially any x86-64 from ~2009+). The install script reads this
# env var, so it applies even though setup-all doesn't pass --march. Override
# with RACCESS_MARCH=native on a homogeneous machine for a faster binary.
export RACCESS_MARCH="${RACCESS_MARCH:-x86-64-v2}"

# Optional: set OFFHOME_CACHE=1 to keep caches/temp off a quota-limited $HOME.
if [ "${OFFHOME_CACHE:-0}" = "1" ]; then
  export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$BASE/.cache}"
  export TMPDIR="${TMPDIR:-$BASE/.tmp}"
  mkdir -p "$XDG_CACHE_HOME" "$TMPDIR"
fi

mkdir -p "$TAUSO_DATA_DIR"
cd "$REPO"
run() { "$MM" run -n "$ENV" "$@"; }

echo ">>> [1/4] setup-all  (genome, bowtie, omics[depmap/mrna/tgcn/attract/riboseq], raccess -march=$RACCESS_MARCH)"
run tauso setup-all -t "$THREADS" --mem-per-thread "$MEM_PER_THREAD"

echo ">>> [2/4] build-cell-context  (cohort, cohort expression, CAI weights, tGCN)"
run tauso build-cell-context

echo ">>> [3/4] assign canonical genes  (process oligo data)"
run python -u -m notebooks.data.OligoAI.assign_canonical_gene

echo ">>> [4/4] index oligo data"
run python -u -m notebooks.utils.data

echo ">>> data setup complete. TAUSO_DATA_DIR=$TAUSO_DATA_DIR"
