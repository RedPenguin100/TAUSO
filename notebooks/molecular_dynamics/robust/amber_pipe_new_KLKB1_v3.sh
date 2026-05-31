#!/bin/bash
# amber_pipe v3 — hardened KLKB1 MD pipeline (drop-in for upstream/).
# Over upstream it adds: resume-from-stage (skip a stage whose restart
# already exists), .done/.failed.<jobid>.<stage> markers, a TERM/INT trap,
# skipping the topology rebuild if .prmtop exists, and --bind-to none on the
# CPU mpirun. The five MD stages share one run_stage helper.
# Usage (unchanged): ./amber_pipe_new_KLKB1_v3.sh <pdb> <n_threads> <nodes>

set -u

if [ $# -ne 3 ]; then
    echo "Usage: $0 <pdb_file> <n_threads> <nodes>"
    exit 1
fi

pdb_file="$1"
n_threads="$2"
nodes="$3"
total_procs=$((n_threads * nodes))

dir_name=$(dirname "$pdb_file")
base_name=$(basename "$pdb_file" .pdb)
output_path="$dir_name/$base_name"

mark_failed() {
    local stage="${1:-unknown}"
    mkdir -p "$output_path" 2>/dev/null
    touch "$output_path/.failed.${SLURM_JOB_ID:-local}.$stage" 2>/dev/null
    echo "[FAIL] base=$base_name stage=$stage jobid=${SLURM_JOB_ID:-local} at $(date -Iseconds)"
}
trap 'mark_failed signal' TERM INT

echo "[START] base=$base_name pdb=$pdb_file jobid=${SLURM_JOB_ID:-local} at $(date -Iseconds)"

# Skip topology rebuild if .prmtop is already there (e.g. from login-node preflight).
if [ ! -s "$output_path/$base_name.prmtop" ]; then
    ./create_files_new_KLKB1.sh "$pdb_file" || { mark_failed create_files; exit 1; }
else
    echo "[SKIP] topology — $base_name.prmtop already exists"
fi

cd "$output_path" || { mark_failed cd; exit 1; }

USE_GPU=${USE_GPU:-0}
if [ "$USE_GPU" -eq 1 ]; then
    amber_exec="pmemd.cuda"
else
    amber_exec="sander"
fi

# Tleap — only run if .prmtop is missing (and create_files just generated _tleap.in).
if [ ! -s "${base_name}.prmtop" ]; then
    echo "[RUN] tleap at $(date -Iseconds)"
    tleap -f "${base_name}_tleap.in" || { mark_failed tleap; exit 1; }
else
    echo "[SKIP] tleap — ${base_name}.prmtop already exists"
fi

echo "Using Amber executable: $amber_exec"

# run_stage <stage> <in_coords> <out_restart> [traj] [ref]
#   Resume-from-stage: skips if <out_restart> already exists.
#   Writes <stage>.nc when [traj] is non-empty and passes -ref when [ref] is set.
#   On failure, drops a .failed marker and exits.
run_stage() {
    local stage=$1 cin=$2 rout=$3 traj=${4:-} ref=${5:-}
    if [ -s "$rout" ]; then
        echo "[SKIP] $stage — $rout already exists"
        return
    fi
    echo "[RUN] $stage at $(date -Iseconds)"
    local a=( -O
        -i "${base_name}_${stage}.in"
        -o "${base_name}_${stage}.out"
        -p "${base_name}.prmtop"
        -c "$cin"
        -r "$rout" )
    [ -n "$traj" ] && a+=( -x "${base_name}_${stage}.nc" )
    [ -n "$ref" ]  && a+=( -ref "$ref" )
    mpirun --bind-to none -np ${total_procs} "$amber_exec".MPI "${a[@]}" \
        || { mark_failed "$stage"; exit 1; }
}

run_stage min1        "${base_name}.rst7"              "${base_name}_min1.ncrst"        ""   "${base_name}.rst7"
run_stage min2        "${base_name}_min1.ncrst"        "${base_name}_min2.ncrst"
run_stage heat        "${base_name}_min2.ncrst"        "${base_name}_heat.ncrst"        traj "${base_name}_min2.ncrst"
run_stage equilibrate "${base_name}_heat.ncrst"        "${base_name}_equilibrate.ncrst" traj "${base_name}_heat.ncrst"
run_stage md          "${base_name}_equilibrate.ncrst" "${base_name}_md.rst7"           traj

# Success — clear stale failure markers and write .done.
rm -f .failed.*
touch .done
echo "[DONE] base=$base_name at $(date -Iseconds)"
