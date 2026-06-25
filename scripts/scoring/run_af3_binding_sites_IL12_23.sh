#!/bin/bash
#SBATCH --job-name=af3_binding_sites_IL12_23
#SBATCH --output=/shared/B3/IL12_23/af3_binding_sites_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --priority=TOP

set -e

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="/shared/B3/IL12_23"
PROTEINS="IL12RB1 IL12p35 IL12p40 IL23p19"
SCRIPTS_ROOT="/home/ubuntu/screening_workflow/scripts"
GENERAL_PYTHON="/home/ubuntu/miniconda3/envs/general/bin/python"

for protein in $PROTEINS; do
    af3_out="${TASK_ROOT}/${protein}/fine_screening/AF3/output"
    af3_dir="${TASK_ROOT}/${protein}/fine_screening/AF3"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing ${protein}..."
    $GENERAL_PYTHON "${SCRIPTS_ROOT}/scoring/af3_scores.py" \
        --af3-results-folder "${af3_out}" \
        --output-dir "${af3_dir}" \
        --binding-sites \
        --distance-threshold 10.0
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done: ${af3_dir}/summary.csv"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All proteins complete."
