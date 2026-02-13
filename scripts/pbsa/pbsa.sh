#!/bin/bash
#SBATCH --job-name=pbsa
#SBATCH --time=480:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yangl@pacagen.com
#SBATCH -C g5.xlarge
#

cd $SLURM_SUBMIT_DIR

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate gmxMMPBSA 

input_parent_dir=$1
range1=$2
range2=$3
output_dir=$4
script_dir=""

for i in {$range1..$range2} ;
do
#input_path=/shared2/B2/stage3_docking/PBSA/MD/target/target_$i/
input_path=$input_parent_dir/$i
if [ -f "$input_path/T298.gro" ]; then
cp -r /home/ubuntu/Applications/gromacs-2023.3/share/top/amber99sb-ildn.ff/ $input_path/.
mkdir $output_dir/$i
cd $output_dir/$i
echo | pwd
env "PATH=$PATH" /home/ubuntu/miniconda3/envs/gmxMMPBSA/bin/gmx_MMPBSA -O -i $script_dir/mmpbsa.in -cs $input_path/T298.tpr -ct $input_path/T298.xtc -ci $input_path/index.ndx -cg 1 13 -cp $input_path/system.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv
fi
done
exit 0
