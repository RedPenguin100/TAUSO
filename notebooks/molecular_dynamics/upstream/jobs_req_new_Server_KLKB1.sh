#!/bin/bash

# Usage: ./jobs_req.sh {dir_path} {n_threads} {nodes} {queue_name}

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

# Iterate over pdb files in the directory
for pdb_file in "$dir_path"/*.pdb; do
    if [ ! -f "$pdb_file" ]; then
        echo "No .pdb files found in $dir_path"
        exit 1
    fi

    dir_name=$(dirname "$pdb_file")
    file_name=$(basename "$pdb_file" .pdb)
    echo "combine"
    output_path="$dir_name/$file_name"
    job_file="${file_name}.sh"

    echo "Creating job file for $file_name..."

    # Create job file (Slurm)
    cat > "$job_file" <<EOF
#!/bin/bash
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${n_threads}
#SBATCH --partition=${queue_name}
#SBATCH -A tamirtul-users_v2
#SBATCH --qos=owner
#SBATCH --mem=250G
#SBATCH --time=2-00:00:00

hostname

# Use your own Miniconda directly (no modules needed)
export PATH="/tamir2/${name}/miniconda3/condabin:\$PATH"
source "/tamir2/${name}/miniconda3/etc/profile.d/conda.sh"

conda activate /tamir2/${name}/miniconda3/envs/amber25

cd /tamir2/${name}/amber
./amber_pipe_new_KLKB1.sh "$output_path" "${n_threads}" "${nodes}"
EOF

    echo "Submitting $job_file..."
    sbatch "$job_file"

    echo "Deleting $job_file..."
    rm "$job_file"

    echo "Job for $file_name submitted and job file removed."
done

echo "All jobs submitted."

