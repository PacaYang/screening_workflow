#!/bin/bash
#SBATCH --job-name=vina_compile
#SBATCH --partition=g24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=14G
#SBATCH --time=02:00:00
#SBATCH --output=/home/yangl_pacagen_com/Projects/B3/IL4_Apr16/logs/vina_compile_%A_%a.out
#SBATCH --error=/home/yangl_pacagen_com/Projects/B3/IL4_Apr16/logs/vina_compile_%A_%a.err
#SBATCH --array=0-3

PROTEINS=(IL4 IL4RA IL13 IL13RA1)
PROTEIN=${PROTEINS[$SLURM_ARRAY_TASK_ID]}
TASK_ROOT=/home/yangl_pacagen_com/Projects/B3/IL4_Apr16

echo "=== vina_compile: ${PROTEIN} (task ${SLURM_ARRAY_TASK_ID}) ==="
echo "Started: $(date)"

cd /home/yangl_pacagen_com/screening_workflow

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

# --- Vina (boltz_test has tqdm/pandas/numpy; vina_test lacks tqdm) ---
echo "--- Vina scoring ---"
conda activate boltz_test
python scripts/scoring/vina_scores.py \
  --vina-results-folder ${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/output \
  --input-dir ${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/input \
  --output-dir ${TASK_ROOT}/${PROTEIN}/fine_screening/Vina

echo "Vina exit code: $?"

# --- Compile ---
echo "--- Compiling scores ---"
conda activate HMSA_test
python scripts/scoring/compile_scores.py \
  --target-dir ${TASK_ROOT}/${PROTEIN} \
  --target-name ${PROTEIN} \
  --output ${TASK_ROOT}/${PROTEIN}/fine_screening/compiled_scores.csv \
  --rfaa-path ${TASK_ROOT}/${PROTEIN}/fine_screening/RoseTTAFold/summary.csv

echo "Compile exit code: $?"
echo "=== Done: $(date) ==="
