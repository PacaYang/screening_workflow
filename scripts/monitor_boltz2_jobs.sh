#!/bin/bash
#
# Monitor Boltz2 Batch Jobs
# This script helps monitor the progress of submitted Boltz2 batch jobs
#

set -e

TASK_ROOT="/home/ubuntu/snake_test"

# Get protein list
get_proteins() {
    echo "JAK1JH1"  # Replace with your protein list
}

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# ============================================================================
# Main monitoring
# ============================================================================

echo "========================================"
echo "Boltz2 Batch Job Monitoring"
echo "========================================"
echo ""

for PROTEIN in $(get_proteins); do
    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/output"
    TOKEN_DIR="${OUTPUT_DIR}/token"
    JOB_LIST="${OUTPUT_DIR}/submitted_jobs.txt"

    if [ ! -f "$JOB_LIST" ]; then
        echo "No jobs submitted yet for protein: ${PROTEIN}"
        continue
    fi

    echo "Protein: ${PROTEIN}"
    echo "----------------------------------------"

    # Count total submitted jobs
    TOTAL_JOBS=$(wc -l < "$JOB_LIST")
    echo "Total batches submitted: ${TOTAL_JOBS}"

    # Count completed batches
    if [ -d "$TOKEN_DIR" ]; then
        COMPLETED=$(find "$TOKEN_DIR" -name "batch_*.done" | wc -l)
    else
        COMPLETED=0
    fi
    echo "Completed batches: ${COMPLETED}"

    # Count compressed outputs
    if [ -d "$OUTPUT_DIR" ]; then
        TAR_FILES=$(find "$OUTPUT_DIR" -maxdepth 1 -name "batch_*.tar.gz" | wc -l)
    else
        TAR_FILES=0
    fi
    echo "Compressed outputs: ${TAR_FILES}"

    # Calculate progress
    if [ $TOTAL_JOBS -gt 0 ]; then
        PERCENT=$((COMPLETED * 100 / TOTAL_JOBS))
        echo "Progress: ${PERCENT}%"
    fi

    echo ""

    # Check SLURM job status
    echo "SLURM Job Status:"
    echo "  Running: $(squeue -u $USER --name="boltz2_${PROTEIN}_*" --states=R | tail -n +2 | wc -l)"
    echo "  Pending: $(squeue -u $USER --name="boltz2_${PROTEIN}_*" --states=PD | tail -n +2 | wc -l)"
    echo ""

    # Show recent completions
    if [ -d "$TOKEN_DIR" ] && [ $COMPLETED -gt 0 ]; then
        echo "Recently completed batches (last 5):"
        find "$TOKEN_DIR" -name "batch_*.done" -printf "%T@ %f\n" | \
            sort -rn | head -5 | cut -d' ' -f2 | sed 's/^/  - /'
        echo ""
    fi

    # Check for failed jobs by looking at error logs
    if [ -d "$OUTPUT_DIR" ]; then
        FAILED_JOBS=$(grep -l "failed\|error\|Error\|ERROR" "${OUTPUT_DIR}"/slurm_batch_*_*.err 2>/dev/null | wc -l)
        if [ $FAILED_JOBS -gt 0 ]; then
            echo "⚠ WARNING: ${FAILED_JOBS} batch(es) may have errors. Check error logs."
            echo ""
        fi
    fi

    echo "========================================"
    echo ""
done

# Overall SLURM queue status
echo "Overall SLURM Queue Status:"
echo "----------------------------------------"
squeue -u $USER -o "%.18i %.9P %.50j %.8T %.10M %.6D %R" | grep -E "JOBID|boltz2_" || echo "No Boltz2 jobs in queue"
echo ""

echo "========================================"
echo "Monitoring completed at: $(date)"
echo "========================================"
