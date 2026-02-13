#!/bin/bash
#
# Monitor DiffDock Docking Jobs
# This script helps monitor the progress of submitted DiffDock docking jobs
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
echo "DiffDock Docking Job Monitoring"
echo "========================================"
echo ""

for PROTEIN in $(get_proteins); do
    INPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/input"
    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/DiffDock/output"
    JOB_LIST="${OUTPUT_DIR}/submitted_jobs.txt"

    if [ ! -f "$JOB_LIST" ]; then
        echo "No jobs submitted yet for protein: ${PROTEIN}"
        continue
    fi

    echo "Protein: ${PROTEIN}"
    echo "----------------------------------------"

    # Count total submitted jobs
    TOTAL_JOBS=$(wc -l < "$JOB_LIST")
    echo "Total parts submitted: ${TOTAL_JOBS}"

    # Count completed parts (by checking .done tokens)
    if [ -d "$OUTPUT_DIR" ]; then
        COMPLETED=$(find "$OUTPUT_DIR" -maxdepth 1 -name "*.done" | wc -l)
    else
        COMPLETED=0
    fi
    echo "Completed parts: ${COMPLETED}"

    # Count total input parts
    if [ -d "$INPUT_DIR" ]; then
        TOTAL_PARTS=$(ls "$INPUT_DIR"/*.csv 2>/dev/null | wc -l)
        echo "Total parts expected: ${TOTAL_PARTS}"
    fi

    # Calculate progress
    if [ $TOTAL_JOBS -gt 0 ]; then
        PERCENT=$((COMPLETED * 100 / TOTAL_JOBS))
        echo "Progress: ${PERCENT}%"
    fi

    echo ""

    # Check SLURM job status
    echo "SLURM Job Status:"
    echo "  Running: $(squeue -u $USER --name="diffdock_${PROTEIN}_*" --states=R 2>/dev/null | tail -n +2 | wc -l)"
    echo "  Pending: $(squeue -u $USER --name="diffdock_${PROTEIN}_*" --states=PD 2>/dev/null | tail -n +2 | wc -l)"
    echo ""

    # Show recently completed parts
    if [ -d "$OUTPUT_DIR" ] && [ $COMPLETED -gt 0 ]; then
        echo "Recently completed parts (last 5):"
        find "$OUTPUT_DIR" -maxdepth 1 -name "*.done" -printf "%T@ %f\n" 2>/dev/null | \
            sort -rn | head -5 | cut -d' ' -f2 | sed 's/^/  - /'
        echo ""
    fi

    # Check for failed jobs by looking at error logs
    if [ -d "$OUTPUT_DIR" ]; then
        # Look for actual failures (not just warnings)
        FAILED_JOBS=0
        for err_file in "${OUTPUT_DIR}"/slurm_part_*_*.err; do
            if [ -f "$err_file" ]; then
                if grep -q "failed\|Failed\|FAILED\|exit code [1-9]" "$err_file" 2>/dev/null; then
                    FAILED_JOBS=$((FAILED_JOBS + 1))
                fi
            fi
        done

        if [ $FAILED_JOBS -gt 0 ]; then
            echo "⚠ WARNING: ${FAILED_JOBS} part(s) may have errors. Check error logs."
            echo ""
        fi
    fi

    # Show sample of output structure
    if [ $COMPLETED -gt 0 ]; then
        echo "Sample output structure:"
        # DiffDock creates directories for each protein_ligand pair
        SAMPLE_DIRS=$(find "$OUTPUT_DIR" -maxdepth 1 -type d -name "${PROTEIN}_*" 2>/dev/null | head -3)
        if [ -n "$SAMPLE_DIRS" ]; then
            echo "  Found $(find "$OUTPUT_DIR" -maxdepth 1 -type d -name "${PROTEIN}_*" 2>/dev/null | wc -l) docking results"
            echo "  Sample directories:"
            echo "$SAMPLE_DIRS" | head -3 | while read dir; do
                echo "    - $(basename "$dir")"
            done
        fi
        echo ""
    fi

    # Estimate total output size
    if [ -d "$OUTPUT_DIR" ]; then
        TOTAL_SIZE=$(du -sh "$OUTPUT_DIR" 2>/dev/null | cut -f1)
        echo "Total output size: ${TOTAL_SIZE}"
        echo ""
    fi

    echo "========================================"
    echo ""
done

# Overall SLURM queue status
echo "Overall SLURM Queue Status:"
echo "----------------------------------------"
squeue -u $USER -o "%.18i %.9P %.50j %.8T %.10M %.6D %R" | grep -E "JOBID|diffdock_" || echo "No DiffDock jobs in queue"
echo ""

echo "========================================"
echo "Monitoring completed at: $(date)"
echo "========================================"
