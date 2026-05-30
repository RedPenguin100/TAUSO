#!/bin/bash
# amber_pipe v3. See ../IMPROVEMENTS.md for what differs from upstream.
# Adds: --bind-to none (§0), .done/.failed markers (§1), resume-from-stage
# (§2), TERM/INT trap (§6), skip if .prmtop exists (§9).
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

# §9: skip topology rebuild if .prmtop is already there (e.g. from login-node preflight)
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

# Tleap — only run if .prmtop is missing (and create_files just generated _tleap.in)
if [ ! -s "${base_name}.prmtop" ]; then
    echo "[RUN] tleap at $(date -Iseconds)"
    tleap -f "${base_name}_tleap.in" || { mark_failed tleap; exit 1; }
else
    echo "[SKIP] tleap — ${base_name}.prmtop already exists"
fi

echo "Using Amber executable: $amber_exec"

# §2: resume-from-stage — skip any stage whose restart output is already present.

# --- Minimization 1 ---
if [ ! -s "${base_name}_min1.ncrst" ]; then
    echo "[RUN] min1 at $(date -Iseconds)"
    mpirun --bind-to none -np ${total_procs} "$amber_exec".MPI -i "${base_name}_min1.in" \
        -o "${base_name}_min1.out" \
        -p "${base_name}.prmtop" \
        -c "${base_name}.rst7" \
        -r "${base_name}_min1.ncrst" \
        -ref "${base_name}.rst7"
    [ $? -ne 0 ] && { mark_failed min1; exit 1; }
else
    echo "[SKIP] min1 — ${base_name}_min1.ncrst already exists"
fi

# --- Minimization 2 ---
if [ ! -s "${base_name}_min2.ncrst" ]; then
    echo "[RUN] min2 at $(date -Iseconds)"
    mpirun --bind-to none -np ${total_procs} "$amber_exec".MPI -i "${base_name}_min2.in" \
        -o "${base_name}_min2.out" \
        -p "${base_name}.prmtop" \
        -c "${base_name}_min1.ncrst" \
        -r "${base_name}_min2.ncrst"
    [ $? -ne 0 ] && { mark_failed min2; exit 1; }
else
    echo "[SKIP] min2 — ${base_name}_min2.ncrst already exists"
fi

# --- Heat ---
if [ ! -s "${base_name}_heat.ncrst" ]; then
    echo "[RUN] heat at $(date -Iseconds)"
    mpirun --bind-to none -np ${total_procs} "$amber_exec".MPI -O -i "${base_name}_heat.in" \
        -o "${base_name}_heat.out" \
        -p "${base_name}.prmtop" \
        -c "${base_name}_min2.ncrst" \
        -r "${base_name}_heat.ncrst" \
        -x "${base_name}_heat.nc" \
        -ref "${base_name}_min2.ncrst"
    [ $? -ne 0 ] && { mark_failed heat; exit 1; }
else
    echo "[SKIP] heat — ${base_name}_heat.ncrst already exists"
fi

# --- Equilibrate ---
if [ ! -s "${base_name}_equilibrate.ncrst" ]; then
    echo "[RUN] equilibrate at $(date -Iseconds)"
    mpirun --bind-to none -np ${total_procs} "$amber_exec".MPI -O -i "${base_name}_equilibrate.in" \
        -o "${base_name}_equilibrate.out" \
        -p "${base_name}.prmtop" \
        -c "${base_name}_heat.ncrst" \
        -r "${base_name}_equilibrate.ncrst" \
        -x "${base_name}_equilibrate.nc" \
        -ref "${base_name}_heat.ncrst"
    [ $? -ne 0 ] && { mark_failed equilibrate; exit 1; }
else
    echo "[SKIP] equilibrate — ${base_name}_equilibrate.ncrst already exists"
fi

# --- Production MD ---
if [ ! -s "${base_name}_md.rst7" ]; then
    echo "[RUN] md at $(date -Iseconds)"
    mpirun --bind-to none -np ${total_procs} "$amber_exec".MPI -O -i "${base_name}_md.in" \
        -o "${base_name}_md.out" \
        -p "${base_name}.prmtop" \
        -c "${base_name}_equilibrate.ncrst" \
        -r "${base_name}_md.rst7" \
        -x "${base_name}_md.nc"
    [ $? -ne 0 ] && { mark_failed md; exit 1; }
else
    echo "[SKIP] md — ${base_name}_md.rst7 already exists"
fi

# §1: success — clear stale failure markers and write .done
rm -f .failed.*
touch .done
echo "[DONE] base=$base_name at $(date -Iseconds)"
