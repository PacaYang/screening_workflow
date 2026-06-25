#!/bin/bash
#SBATCH --job-name=peptide_af3_validate
#SBATCH --output=/home/yangl_pacagen_com/peptide_test/slurm_%j.out
#SBATCH --error=/home/yangl_pacagen_com/peptide_test/slurm_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=g24

set -euo pipefail

OUTPUT_DIR=/home/yangl_pacagen_com/peptide_test
PROTEIN_NAME=2O0O

echo "Starting AF3 validation pipeline at $(date)"

# Step 4: Submit AF3 jobs
echo "Step 4: Submitting AF3 jobs..."
bash /home/yangl_pacagen_com/screening_workflow/scripts/peptide_generation/step4_submit_af3.sh \
  --output-dir $OUTPUT_DIR \
  --protein-name $PROTEIN_NAME \
  --n-batches 2

# Wait for AF3 jobs to complete
echo "Waiting for AF3 jobs..."
while squeue -u $USER -n "${PROTEIN_NAME}_af3_*" --noheader 2>/dev/null | grep -q .; do
    sleep 60
done

# Step 5: Collect results
echo "Step 5: Collecting results..."
bash /home/yangl_pacagen_com/screening_workflow/scripts/peptide_generation/step5_collect_streaming.sh \
  --output-dir $OUTPUT_DIR \
  --protein-name $PROTEIN_NAME \
  --max-iterations 10

# Step 6: Finalize ranking
echo "Step 6: Finalizing ranking..."
python3 /home/yangl_pacagen_com/screening_workflow/scripts/peptide_generation/step6_finalize_ranking.py \
  --output-dir $OUTPUT_DIR \
  --protein-name $PROTEIN_NAME

echo "Pipeline complete at $(date)"
echo "Results: ${OUTPUT_DIR}/05_results/ranked_peptides.csv"
