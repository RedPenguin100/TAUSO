#!/bin/bash

# Usage: ./jobs_req.sh {name} {dir_path}

if [ $# -ne 3 ]; then
    echo "Usage: $0 {dir_path} {n_threads} {nodes}"
    exit 1
fi
dir_path="$1"
n_threads="$2"
nodes="$3"

name="$(whoami)"
echo ${name}


# Iterate over pdb files in the directory
for pdb_file in "$dir_path"/*.pdb; do
    if [ ! -f "$pdb_file" ]; then
        echo "No .pdb files found in $dir_path"
        exit 1
    fi

    # Get directory path of the input file
    dir_name=$(dirname "$pdb_file")

    # Extract the base name without extension
    file_name=$(basename "$pdb_file" .pdb)

    # Combine both to stay inside the same path
    output_path="$dir_name/$file_name"
    job_file="${file_name}.sh"

    echo "Creating job file for $file_name..."

    # Create job file
    cat > "$job_file" <<EOF
#!/bin/bash
#PBS -l select=${nodes}:ncpus=${n_threads}:mpiprocs=${n_threads}

hostname

module load miniconda/miniconda3-2023-environmentally
conda activate /tamir2/${name}/miniconda3/envs/amber25

cd /tamir2/${name}/amber
./amber_pipe.sh "$output_path" "${n_threads}" "${nodes}"
EOF

    # Submit the job
    echo "Submitting $job_file..."
    qsub -q tamirQ "$job_file"

    # Delete job file
    echo "Deleting $job_file..."
    rm "$job_file"

    echo "Job for $file_name submitted and job file removed."
done

echo "All jobs submitted."

