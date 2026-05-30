#!/bin/bash
# jobs_req — robust v3 (public partition).
#
# Improvements over upstream `jobs_req_new_Server_KLKB1_P.sh`:
#   §3  --mem dropped from 250G to 16G. The original ask was 50–250× the
#       actual sander.MPI footprint for a 20-mer in OPC. Reducing it lets
#       jobs land on smaller nodes and bumps queue throughput. If a real
#       benchmark shows MaxRSS > 12 GB, raise this.
#   §5  Tleap preflight runs on the login node before sbatch. Malformed
#       PDBs are logged to <batch>/bad_pdbs.txt and skipped — they no
#       longer eat a SLURM allocation. Side effect: the SLURM job sees
#       .prmtop already present and skips tleap (handled by v3 amber_pipe §9).
#   §6  --requeue: SLURM auto-resubmits on node failure / preemption.
#       Pair with v3 amber_pipe's resume-from-stage and signal trap.
#
# Invokes amber_pipe_new_KLKB1_v3.sh and create_files_new_KLKB1.sh; both
# must sit next to this script in the same directory (i.e. /tamir2/$USER/amber).
#
# Conda env is hardcoded to Ariella's amber25 because she owns the
# canonical install; replace if you have your own.
#
# Usage (unchanged from upstream): ./jobs_req_new_Server_KLKB1_P_v3.sh \
#       <batch_dir> <n_threads> <nodes> <partition>

if [ $# -ne 4 ]; then
    echo "Usage: $0 {dir_path} {n_threads} {nodes} {queue_name}"
    exit 1
fi
dir_path="$1"
n_threads="$2"
nodes="$3"
queue_name="$4"

name="$(whoami)"
echo "${name}"
echo "start"

# Conda env (Ariella's amber25)
CONDA_ROOT="/tamir2/nouman/miniconda3"
CONDA_ENV="$CONDA_ROOT/envs/amber25"

# §5 preflight bad-PDB log
bad_pdbs_log="$dir_path/bad_pdbs.txt"

for pdb_file in "$dir_path"/*.pdb; do
    if [ ! -f "$pdb_file" ]; then
        echo "No .pdb files found in $dir_path"
        exit 1
    fi

    dir_name=$(dirname "$pdb_file")
    file_name=$(basename "$pdb_file" .pdb)
    output_path="$dir_name/$file_name"
    job_file="${file_name}.sh"

    echo "[PREFLIGHT] $file_name"

    # Make .in files (idempotent — create_files mkdir -p's and overwrites .in files).
    ./create_files_new_KLKB1.sh "$pdb_file" > /dev/null 2>&1

    # §5 tleap on login — only if .prmtop is missing (§9-aware: skips on retries).
    if [ ! -s "$output_path/${file_name}.prmtop" ]; then
        (
            source "$CONDA_ROOT/etc/profile.d/conda.sh"
            conda activate "$CONDA_ENV"
            cd "$output_path"
            tleap -f "${file_name}_tleap.in"
        ) > "$output_path/tleap_preflight.log" 2>&1
    fi

    if [ ! -s "$output_path/${file_name}.prmtop" ]; then
        echo "[BAD] $file_name — tleap failed on login; logged + skipped"
        echo "$(date -Iseconds)	$pdb_file" >> "$bad_pdbs_log"
        continue
    fi

    echo "Creating job file for $file_name..."

    # Create SLURM job file
    cat > "$job_file" <<EOF
#!/bin/bash
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${n_threads}
#SBATCH --partition=${queue_name}
#SBATCH -A tamirtul-users_v2
#SBATCH --qos=public
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --requeue

hostname
date -Iseconds

export PATH="$CONDA_ROOT/condabin:\$PATH"
source "$CONDA_ROOT/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

cd /tamir2/${name}/amber
./amber_pipe_new_KLKB1_v3.sh "$output_path" "${n_threads}" "${nodes}"
EOF

    echo "Submitting $job_file..."
    sbatch "$job_file"

    echo "Deleting $job_file..."
    rm "$job_file"

    echo "Job for $file_name submitted and job file removed."
done

echo "All jobs submitted."
