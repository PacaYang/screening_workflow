#!/bin/bash 

# Check the forcefield type: amber99sb-ildn.ff
# include the ligand itp, gaff_atomtype.itp in the system_EM

sdf=$1
output_dir=$2
input_gro=$3
input_EM_top=$4
script_dir=$5

python=/home/ubuntu/miniconda3/envs/gmxMMPBSA/bin/python
obabel=/home/ubuntu/miniconda3/envs/gmxMMPBSA/bin/obabel
acpype=/home/ubuntu/miniconda3/envs/gmxMMPBSA/bin/acpype

# change to working directory
if [ -d $output_dir ]; then
        :
else
        mkdir -p $output_dir
fi

cd $output_dir

if [ ! -f T298.gro ]; then

$python $script_dir/write_mol.py -input $sdf
$obabel ligand.mol -O ligand.mol2 -h
$acpype -i ligand.mol2 -b ligand

echo "================================ Ligand ff generated =========================================="

# Count atoms
n_prot=$(awk 'NR==2' $input_gro)
n_lig=$(awk 'NR==2' ligand.acpype/ligand_GMX.gro)
total_atoms=$((n_prot + n_lig))

# Extract coordinate lines
head -n -1 $input_gro | tail -n +3 > prot_coords.txt
head -n -1 ligand.acpype/ligand_GMX.gro | tail -n +3 > lig_coords.txt

# Extract box vector (last line)
tail -n 1 $input_gro > box.txt

# Combine
(echo "Protein-Ligand complex"
 echo "$total_atoms"
 cat prot_coords.txt
 cat lig_coords.txt
 cat box.txt) > complex.gro

awk '/\[ *atomtypes *\]/ {p=1; print; next} /^\[/ && p {exit} p' ligand.acpype/ligand_GMX.itp > gaff_atomtypes.itp
awk 'BEGIN{p=1} /\[ *atomtypes *\]/ {p=0; next} /^\[/ && !p {p=1} p' ligand.acpype/ligand_GMX.itp > ligand.itp


source /home/ubuntu/Applications/gromacs-2025.3/bin/GMXRC
gmx=/home/ubuntu/Applications/gromacs-2025.3/bin/gmx

$gmx editconf -f complex.gro -o complex.gro -d 2.5
cp $input_EM_top system.top

$gmx solvate -cp complex.gro -cs spc216 -p system -o sol
$gmx grompp -f $script_dir/MDP/em -c sol -o ions -p system
echo -e "SOL\n" | gmx genion -s ions -neutral -p system -conc 0.15 -o ions

echo -e "r UNL\nq\n" | $gmx make_ndx -f ions.gro

$gmx grompp -f $script_dir/MDP/em -p system -c ions.gro -o em
$gmx mdrun -deffnm em -v

$gmx grompp -f $script_dir/MDP/md1 -p system -c em.gro -o md1
$gmx mdrun -deffnm md1 -v

$gmx grompp -f $script_dir/MDP/md2 -p system -c md1.gro -o md2 -maxwarn 1
$gmx mdrun -deffnm md2 -v 

$gmx grompp -f $script_dir/MDP/T298 -p system -c md2.gro -o T298 -maxwarn 1
$gmx mdrun -deffnm T298 -cpi T298 -v

fi
