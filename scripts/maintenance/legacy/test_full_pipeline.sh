#!/bin/bash
#
# Test script for RoseTTAFold and Vina on IL6RA
# Uses run_full_pipeline.sh with skip flags to isolate these two methods.
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TASK_ROOT="/home/yangl_pacagen_com/Projects/B3/IL6"
PROTEIN="IL6RA"

echo "============================================"
echo "  Test: RoseTTAFold + Vina for ${PROTEIN}"
echo "============================================"
echo ""

# ── Step 1: Clean up cached files from previous failed Vina run ──
echo "Cleaning up cached Vina files from previous failed run..."
VINA_OUT="${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/output/input_part_0"
rm -f "${VINA_OUT}/IL6RAFH.pdb"
rm -f "${VINA_OUT}/IL6RAFH.txt"
rm -f "${VINA_OUT}/IL6RA_minimized.pdb"
rm -f "${VINA_OUT}/IL6RAFH_box"*
rm -f "${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/output/input_part_0.done"
echo "  Done."

# ── Step 2: Clean up failed RoseTTAFold prefold state ──
echo "Cleaning up failed RoseTTAFold prefold state..."
RFAA_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/RoseTTAFold"
rm -f "${RFAA_DIR}/submitted_jobs.txt"
# Don't remove MSA outputs (hhblits, a3m) — they're expensive to regenerate
# Only remove the fold token so prefold re-runs
rm -f "${RFAA_DIR}/protein_folding/output/protein_fold.done"
echo "  Done."
echo ""

# ── Step 3: Run the pipeline (stages 4+5 only, RoseTTAFold + Vina only) ──
echo "Launching run_full_pipeline.sh..."
echo ""

bash "${SCRIPT_DIR}/run_full_pipeline.sh" \
    --task-root "${TASK_ROOT}" \
    --proteins "${PROTEIN}" \
    --start-from 4 \
    --stop-after 5 \
    --skip-graphdta \
    --skip-hmsa \
    --skip-colddta \
    --skip-druglamp \
    --skip-conplex \
    --skip-af3 \
    --skip-boltz2 \
    --skip-diffdock \
    --skip-md-pbsa \
    --poll-interval 60
