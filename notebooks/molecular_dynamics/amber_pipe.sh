#!/bin/bash
# Check if a PDB file was provided

if [ $# -ne 3 ]; then
    echo "Usage: $0 <pdb_file> <n_threads> <nodes>"
    exit 1
fi

pdb_file="$1"
n_threads="$2"
nodes="$3"
total_procs=$((n_threads * nodes))


echo "Using PDB file: $pdb_file"


# Get directory path of the input file
dir_name=$(dirname "$pdb_file")

# Extract the base name without extension
base_name=$(basename "$pdb_file" .pdb)

# Combine both to stay inside the same path
output_path="$dir_name/$base_name"

./create_files.sh "$pdb_file"

cd $output_path

# Example: check for GPU or CPU
# Set this manually or detect automatically (simple example: check for pmemd.cuda)
USE_GPU=${USE_GPU:-0}  # default to 0 (CPU)

if [ "$USE_GPU" -eq 1 ]; then
    amber_exec="pmemd.cuda"
else
    amber_exec="sander"
fi



tleap -f "${base_name}_tleap.in"

echo "Using Amber executable: $amber_exec"

# --- Minimization 1 ---
mpirun -np ${total_procs} "$amber_exec".MPI -i "${base_name}_min1.in" \
    -o "${base_name}_min1.out" \
    -p "${base_name}.prmtop" \
    -c "${base_name}.rst7" \
    -r "${base_name}_min1.ncrst" \
    -ref "${base_name}.rst7"

if [ $? -ne 0 ]; then
    echo "Error in minimization 1, exiting"
    exit 1
fi

# --- Minimization 2 ---
mpirun -np ${total_procs} "$amber_exec".MPI -i "${base_name}_min2.in" \
    -o "${base_name}_min2.out" \
    -p "${base_name}.prmtop" \
    -c "${base_name}_min1.ncrst" \
    -r "${base_name}_min2.ncrst"

if [ $? -ne 0 ]; then
    echo "Error in minimization 2, exiting"
    exit 1
fi

# --- MD 1 ---
mpirun -np ${total_procs} "$amber_exec".MPI -O -i "${base_name}_md1.in" \
    -o "${base_name}_md1.out" \
    -p "${base_name}.prmtop" \
    -c "${base_name}_min2.ncrst" \
    -r "${base_name}_md1.ncrst" \
    -x "${base_name}_md1.nc" \
    -ref "${base_name}_min2.ncrst"

if [ $? -ne 0 ]; then
    echo "Error in MD1, exiting"
    exit 1
fi

# --- MD 2 ---
mpirun -np ${total_procs} "$amber_exec".MPI -O -i "${base_name}_md2.in" \
    -o "${base_name}_md2.out" \
    -p "${base_name}.prmtop" \
    -c "${base_name}_md1.ncrst" \
    -r "${base_name}_md2.ncrst" \
    -x "${base_name}_md2.nc"

if [ $? -ne 0 ]; then
    echo "Error in MD2, exiting"
    exit 1
fi
