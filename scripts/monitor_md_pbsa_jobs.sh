#!/bin/bash
#
# Monitor MD + PBSA Combined Jobs
# This script helps monitor the progress of submitted MD+PBSA jobs
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
echo "MD + PBSA Job Monitoring"
echo "========================================"
echo ""

for PROTEIN in $(get_proteins); do
    PBSA_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/PBSA"
    MD_DIR="${PBSA_DIR}/MD"
    PBSA_RESULTS_DIR="${PBSA_DIR}/PBSA"
    JOB_LIST="${PBSA_DIR}/submitted_jobs.txt"

    if [ ! -f "$JOB_LIST" ]; then
        echo "No jobs submitted yet for protein: ${PROTEIN}"
        continue
    fi

    echo "Protein: ${PROTEIN}"
    echo "----------------------------------------"

    # Count total submitted jobs
    TOTAL_JOBS=$(wc -l < "$JOB_LIST")
    echo "Total jobs submitted: ${TOTAL_JOBS}"

    # Count completed MD steps
    if [ -d "$MD_DIR" ]; then
        MD_COMPLETED=$(find "$MD_DIR" -mindepth 2 -maxdepth 2 -name "token.done" | wc -l)
    else
        MD_COMPLETED=0
    fi
    echo "MD completed: ${MD_COMPLETED}"

    # Count completed PBSA steps
    if [ -d "$PBSA_RESULTS_DIR" ]; then
        PBSA_COMPLETED=$(find "$PBSA_RESULTS_DIR" -mindepth 2 -maxdepth 2 -name "token.done" | wc -l)
    else
        PBSA_COMPLETED=0
    fi
    echo "PBSA completed: ${PBSA_COMPLETED}"

    # Calculate progress
    if [ $TOTAL_JOBS -gt 0 ]; then
        MD_PERCENT=$((MD_COMPLETED * 100 / TOTAL_JOBS))
        PBSA_PERCENT=$((PBSA_COMPLETED * 100 / TOTAL_JOBS))
        echo "MD Progress: ${MD_PERCENT}%"
        echo "PBSA Progress: ${PBSA_PERCENT}%"
    fi

    # Jobs in MD phase (MD done but PBSA not done)
    IN_PBSA_PHASE=$((MD_COMPLETED - PBSA_COMPLETED))
    if [ $IN_PBSA_PHASE -gt 0 ]; then
        echo "In PBSA phase: ${IN_PBSA_PHASE}"
    fi

    echo ""

    # Check SLURM job status
    echo "SLURM Job Status:"
    RUNNING=$(squeue -u $USER --name="md_pbsa_${PROTEIN}_*" --states=R 2>/dev/null | tail -n +2 | wc -l)
    PENDING=$(squeue -u $USER --name="md_pbsa_${PROTEIN}_*" --states=PD 2>/dev/null | tail -n +2 | wc -l)
    echo "  Running: ${RUNNING}"
    echo "  Pending: ${PENDING}"
    echo ""

    # Show recently completed MD steps
    if [ -d "$MD_DIR" ] && [ $MD_COMPLETED -gt 0 ]; then
        echo "Recently completed MD (last 5):"
        find "$MD_DIR" -mindepth 2 -maxdepth 2 -name "token.done" -printf "%T@ %h\n" 2>/dev/null | \
            sort -rn | head -5 | cut -d' ' -f2 | xargs -I {} basename {} | sed 's/^/  - /'
        echo ""
    fi

    # Show recently completed PBSA steps
    if [ -d "$PBSA_RESULTS_DIR" ] && [ $PBSA_COMPLETED -gt 0 ]; then
        echo "Recently completed PBSA (last 5):"
        find "$PBSA_RESULTS_DIR" -mindepth 2 -maxdepth 2 -name "token.done" -printf "%T@ %h\n" 2>/dev/null | \
            sort -rn | head -5 | cut -d' ' -f2 | xargs -I {} basename {} | sed 's/^/  - /'
        echo ""
    fi

    # Check for failed jobs by looking at error logs
    if [ -d "$PBSA_DIR" ]; then
        FAILED_JOBS=0
        SKIPPED_JOBS=0
        for out_file in "${PBSA_DIR}"/slurm_${PROTEIN}_*_*.out; do
            if [ -f "$out_file" ]; then
                # Check for skipped jobs (missing SDF)
                if grep -q "SKIPPING JOB: SDF file not found" "$out_file" 2>/dev/null; then
                    SKIPPED_JOBS=$((SKIPPED_JOBS + 1))
                fi
            fi
        done

        for err_file in "${PBSA_DIR}"/slurm_${PROTEIN}_*_*.err; do
            if [ -f "$err_file" ]; then
                if grep -q "failed\|Failed\|FAILED\|exit code [1-9]" "$err_file" 2>/dev/null; then
                    FAILED_JOBS=$((FAILED_JOBS + 1))
                fi
            fi
        done

        if [ $SKIPPED_JOBS -gt 0 ]; then
            echo "ℹ INFO: ${SKIPPED_JOBS} job(s) skipped (missing DiffDock SDF files)"
        fi

        if [ $FAILED_JOBS -gt 0 ]; then
            echo "⚠ WARNING: ${FAILED_JOBS} job(s) may have errors. Check error logs."
        fi

        if [ $SKIPPED_JOBS -gt 0 ] || [ $FAILED_JOBS -gt 0 ]; then
            echo ""
        fi
    fi

    # Show sample MD output
    if [ $MD_COMPLETED -gt 0 ]; then
        echo "Sample MD output:"
        SAMPLE_MD=$(find "$MD_DIR" -mindepth 1 -maxdepth 1 -type d | head -1)
        if [ -n "$SAMPLE_MD" ]; then
            echo "  Directory: $(basename "$SAMPLE_MD")"
            if [ -f "${SAMPLE_MD}/T298.xtc" ]; then
                echo "  Trajectory: $(du -sh "${SAMPLE_MD}/T298.xtc" 2>/dev/null | cut -f1)"
            fi
            echo "  Files: $(ls "$SAMPLE_MD" 2>/dev/null | wc -l)"
        fi
        echo ""
    fi

    # Show sample PBSA output
    if [ $PBSA_COMPLETED -gt 0 ]; then
        echo "Sample PBSA output:"
        SAMPLE_PBSA=$(find "$PBSA_RESULTS_DIR" -mindepth 1 -maxdepth 1 -type d | head -1)
        if [ -n "$SAMPLE_PBSA" ]; then
            echo "  Directory: $(basename "$SAMPLE_PBSA")"
            if [ -f "${SAMPLE_PBSA}/FINAL_RESULTS_MMPBSA.dat" ]; then
                echo "  Results: FINAL_RESULTS_MMPBSA.dat (present)"
            fi
            if [ -f "${SAMPLE_PBSA}/FINAL_RESULTS_MMPBSA.csv" ]; then
                echo "  CSV: FINAL_RESULTS_MMPBSA.csv (present)"
            fi
        fi
        echo ""
    fi

    # Estimate total output size
    if [ -d "$PBSA_DIR" ]; then
        TOTAL_SIZE=$(du -sh "$PBSA_DIR" 2>/dev/null | cut -f1)
        echo "Total output size: ${TOTAL_SIZE}"
        echo ""
    fi

    # Show runtime statistics for completed jobs
    if [ $PBSA_COMPLETED -gt 0 ]; then
        echo "Runtime statistics (from logs):"
        COMPLETED_LOGS=($(ls -t "${PBSA_DIR}"/slurm_${PROTEIN}_*_*.out 2>/dev/null | head -5))
        if [ ${#COMPLETED_LOGS[@]} -gt 0 ]; then
            for log_file in "${COMPLETED_LOGS[@]}"; do
                if grep -q "Job completed at:" "$log_file" 2>/dev/null; then
                    JOB_NAME=$(basename "$log_file" | sed 's/slurm_//' | sed 's/_[0-9]*\.out//')
                    START_TIME=$(grep "Job started at:" "$log_file" 2>/dev/null | sed 's/.*at: //')
                    END_TIME=$(grep "Job completed at:" "$log_file" 2>/dev/null | sed 's/.*at: //')
                    if [ -n "$START_TIME" ] && [ -n "$END_TIME" ]; then
                        echo "  ${JOB_NAME}: ${START_TIME} → ${END_TIME}"
                    fi
                fi
            done
        fi
        echo ""
    fi

    echo "========================================"
    echo ""
done

# Overall SLURM queue status
echo "Overall SLURM Queue Status:"
echo "----------------------------------------"
squeue -u $USER -o "%.18i %.9P %.50j %.8T %.10M %.6D %R" | grep -E "JOBID|md_pbsa_" || echo "No MD+PBSA jobs in queue"
echo ""

echo "========================================"
echo "Monitoring completed at: $(date)"
echo "========================================"
