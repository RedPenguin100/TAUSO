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
loadamberparams /tamir2/nouman/amber/modXNA-main/dat/frcmod.modxna
loadoff /tamir2/nouman/amber/modXNA-main/PSA.lib
loadoff /tamir2/nouman/amber/modXNA-main/PSC.lib
loadoff /tamir2/nouman/amber/modXNA-main/PSG.lib
loadoff /tamir2/nouman/amber/modXNA-main/PST.lib
loadoff /tamir2/nouman/amber/modXNA-main/FC5.lib
loadoff /tamir2/nouman/amber/modXNA-main/FMA.lib
loadoff /tamir2/nouman/amber/modXNA-main/FMG.lib
loadoff /tamir2/nouman/amber/modXNA-main/FMT.lib
loadoff /tamir2/nouman/amber/modXNA-main/MC5.lib
loadoff /tamir2/nouman/amber/modXNA-main/PC5.lib
loadoff /tamir2/nouman/amber/modXNA-main/PMA.lib
loadoff /tamir2/nouman/amber/modXNA-main/PMG.lib
loadoff /tamir2/nouman/amber/modXNA-main/PMT.lib
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
RES 1 40
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

echo -e "heating: 25 kcal/mol
restraints on DNA, V=const, 100ps
 &cntrl
imin=0, ntx=1, ntpr=500,
ntwr=500, ntwx=500, ntwe=500,
nscm=5000,
ntf=2, ntc=2,
ntb=1, ntp=0,
nstlim=50000, t=0.0, dt=0.002,
cut=10.0,
tempi=100.0, ntt=1,
ntr=1,nmropt=1,
restraintmask=':1-40',
restraint_wt=25.0,
&end
&wt type='TEMP0',istep1=0,
istep2=5000, value1=100.0,value2=310.0,
&end
&wt type='TEMP0',istep1=5001,
istep2=50000, value1=310.0,value2=310.0,
&end
&wt type='END', &end" > "${base_name}_heat.in"

echo -e "equilibrate: 0.5 kcal/mol
restraints on DNA, 500ps, constant T=310K and P=1ATM coupling = 0.2
 &cntrl
imin=0, ntx=5, ntpr=500,
ntwr=500, ntwx=500,
ntwe=500,
nscm=5000,
ntf=2, ntc=2,
ntb=2, ntp=1, tautp=0.2,
taup=0.2,
nstlim=25000,
t=0.0, dt=0.002,
cut=10.0,
ntt=1,
ntr=1,
irest=1,
tempi = 310.0, temp0 = 310.0,
restraintmask=':1-40',
restraint_wt=0.5,
&end
&ewald
ew_type = 0, skinnb = 1.0,
&end" > "${base_name}_equilibrate.in"


echo -e "Main Production: 2ns MD
 &cntrl
  imin = 0, irest = 1, ntx = 7,
  ntb = 2, pres0 = 1.0, ntp = 1,
  nscm = 1000,
  taup = 2.0,
  cut = 10.0,
  ntc = 2, ntf = 2,
  tempi = 310.0, temp0 = 310.0,
  ntt = 3, gamma_ln = 1.0,
  nstlim = 1000000, dt = 0.002,
  ntpr = 1000, ntwx = 1000, ntwr = 1000, ntwe = 1000,
  iwrap = 1, ioutfm   = 1
 /" > "${base_name}_md.in"



