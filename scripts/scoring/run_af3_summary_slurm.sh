#!/bin/bash
#SBATCH --job-name=af3_summary
#SBATCH --partition=g24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=14G
#SBATCH --time=04:00:00
#SBATCH --output=/home/yangl_pacagen_com/Projects/B3/IL4_Apr16/logs/af3_summary_%A_%a.out
#SBATCH --error=/home/yangl_pacagen_com/Projects/B3/IL4_Apr16/logs/af3_summary_%A_%a.err
#SBATCH --array=0-3

PROTEINS=(IL4 IL4RA IL13 IL13RA1)
PROTEIN=${PROTEINS[$SLURM_ARRAY_TASK_ID]}

TASK_ROOT=/home/yangl_pacagen_com/Projects/B3/IL4_Apr16
AF3_DIR=${TASK_ROOT}/${PROTEIN}/fine_screening/AF3

echo "[$(date)] Starting af3_scores.py for ${PROTEIN}"
echo "Input:  ${AF3_DIR}/output"
echo "Output: ${AF3_DIR}/summary.csv"

cd /home/yangl_pacagen_com/screening_workflow

python scripts/scoring/af3_scores.py \
  --af3-results-folder ${AF3_DIR}/output \
  --output-dir ${AF3_DIR}

echo "[$(date)] Done for ${PROTEIN}"
wc -l ${AF3_DIR}/summary.csv
