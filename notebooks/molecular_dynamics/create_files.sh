#!/bin/bash
# Check if a PDB file was provided
if [ -z "$1" ]; then
    echo "Usage: $0 <pdb_file>"
    exit 1
fi

pdb_file="$1"
echo "Using PDB file: $pdb_file"

# Get directory path of the input file
dir_name=$(dirname "$pdb_file")

# Extract the base name without extension
base_name=$(basename "$pdb_file" .pdb)

# Combine both to stay inside the same path
output_path="$dir_name/$base_name"


# Create a directory with this name if it doesn't exist
echo "$output_path"
mkdir -p "$output_path"

cd "$output_path"

echo -e "source leaprc.DNA.OL21
source leaprc.RNA.OL3
source leaprc.water.opc
dna1 = loadpdb ../$base_name.pdb
addions dna1 Na+ 0
solvateoct dna1 OPCBOX 10.0
addionsrand dna1 Na+ 150
addionsrand dna1 Cl- 150
saveamberparm dna1 ${base_name}.prmtop ${base_name}.rst7
quit" > "${base_name}_tleap.in"

echo -e "polyA-polyT 10-mer: initial minimization solvent + ions
 &cntrl
  imin   = 1,
  maxcyc = 1000,
  ncyc   = 500,
  ntb    = 1,
  ntr    = 1,
  cut    = 10.0
 /
Hold the DNA fixed
500.0
RES 1 20
END
END" > "${base_name}_min1.in"

echo -e "polyA-polyT 10-mer: initial minimization whole system
 &cntrl
  imin   = 1,
  maxcyc = 2500,
  ncyc   = 1000,
  ntb    = 1,
  ntr    = 0,
  cut    = 10.0
 /" > "${base_name}_min2.in"

echo -e "polyA-polyT 10-mer: 20ps MD with res on DNA
 &cntrl
  imin   = 0,
  irest  = 0,
  ntx    = 1,
  ntb    = 1,
  cut    = 10.0,
  ntr    = 1,
  ntc    = 2,
  ntf    = 2,
  tempi  = 0.0,
  temp0  = 310.0,
  ntt    = 3,
  gamma_ln = 1.0,
  nstlim = 10000, dt = 0.002
  ntpr = 100, ntwx = 100, ntwr = 1000
 /
Keep DNA fixed with weak restraints
10.0
RES 1 20
END
END" > "${base_name}_md1.in"

echo -e "polyA-polyT 10-mer: 100ps MD
 &cntrl
  imin = 0, irest = 1, ntx = 7,
  ntb = 2, pres0 = 1.0, ntp = 1,
  taup = 2.0,
  cut = 10.0, ntr = 0,
  ntc = 2, ntf = 2,
  tempi = 310.0, temp0 = 310.0,
  ntt = 3, gamma_ln = 1.0,
  nstlim = 500000, dt = 0.002,
  ntpr = 100, ntwx = 100, ntwr = 1000
 /" > "${base_name}_md2.in"



