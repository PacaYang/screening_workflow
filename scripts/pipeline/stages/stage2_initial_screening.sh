#!/bin/bash
#
# Stage 2: Initial Screening
# Description: Run 5 initial screening methods (GraphDTA, HMSA, ColdDTA, DrugLAMP, ConPLex)
# Inputs: ${TASK_ROOT}/${PROTEIN}/initial_screening/inputs/
# Outputs: ${TASK_ROOT}/${PROTEIN}/initial_screening/${METHOD}/prediction_*.csv
# Dependencies: Stage 1
#

set -e

# Source libraries
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../lib/pipeline_utils.sh"
source "${SCRIPT_DIR}/../lib/job_monitor.sh"
source "${SCRIPT_DIR}/../lib/state_manager.sh"
source "${SCRIPT_DIR}/../lib/logger.sh"

# Default values
TASK_ROOT=""
PROTEINS=""
DRY_RUN=0
POLL_INTERVAL=300

# Skip flags
SKIP_GRAPHDTA=0
SKIP_HMSA=0
SKIP_COLDDTA=0
SKIP_DRUGLAMP=0
SKIP_CONPLEX=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --task-root) TASK_ROOT="$2"; shift 2 ;;
        --proteins) PROTEINS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --poll-interval) POLL_INTERVAL="$2"; shift 2 ;;
        --skip-graphdta) SKIP_GRAPHDTA=1; shift ;;
        --skip-hmsa) SKIP_HMSA=1; shift ;;
        --skip-colddta) SKIP_COLDDTA=1; shift ;;
        --skip-druglamp) SKIP_DRUGLAMP=1; shift ;;
        --skip-conplex) SKIP_CONPLEX=1; shift ;;
        --state-file) STATE_FILE="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Validate required arguments
if [ -z "$TASK_ROOT" ] || [ -z "$PROTEINS" ]; then
    log_error "Missing required arguments"
    exit 1
fi

# Main execution
main() {
    log_stage_start "Stage 2: Initial Screening"

    # Update state to running
    [ -n "$STATE_FILE" ] && update_stage_status "2" "running"

    # Export environment variables
    export_pipeline_env "$TASK_ROOT" "$PROTEINS"

    # Submit all enabled methods
    local methods=()
    [ "$SKIP_GRAPHDTA" -eq 0 ] && methods+=("GraphDTA")
    [ "$SKIP_HMSA" -eq 0 ] && methods+=("HMSA")
    [ "$SKIP_COLDDTA" -eq 0 ] && methods+=("ColdDTA")
    [ "$SKIP_DRUGLAMP" -eq 0 ] && methods+=("DrugLAMP")
    [ "$SKIP_CONPLEX" -eq 0 ] && methods+=("ConPLex")

    log_info "Submitting ${#methods[@]} screening methods: ${methods[*]}"

    # Submit jobs
    [ "$SKIP_GRAPHDTA" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../initial_screening/run_graphdta.sh"
    [ "$SKIP_HMSA" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../initial_screening/run_hmsa.sh"
    [ "$SKIP_COLDDTA" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../initial_screening/run_colddta.sh"
    [ "$SKIP_DRUGLAMP" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../initial_screening/run_druglamp.sh"
    [ "$SKIP_CONPLEX" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../initial_screening/run_conplex.sh"

    if [ "$DRY_RUN" -eq 1 ]; then
        log_info "[DRY RUN] Would wait for SLURM jobs to complete"
        return 0
    fi

    # Collect all job files
    local job_files=()
    for method in "${methods[@]}"; do
        for protein in $PROTEINS; do
            local job_file="${TASK_ROOT}/${protein}/initial_screening/${method}/job_ids.txt"
            if [ -f "$job_file" ]; then
                job_files+=("$job_file")
            fi
        done
    done

    # Wait for all jobs to complete
    if [ ${#job_files[@]} -gt 0 ]; then
        log_info "Waiting for ${#job_files[@]} job files to complete"
        wait_for_jobs_with_progress "$POLL_INTERVAL" "${job_files[@]}"
    fi

    # Verify predictions exist
    local ok=1
    for protein in $PROTEINS; do
        for method in "${methods[@]}"; do
            local pred="${TASK_ROOT}/${protein}/initial_screening/${method}"
            local found
            found=$(find "$pred" -maxdepth 1 -name "prediction_*.csv" 2>/dev/null | wc -l)
            if [ "$found" -gt 0 ]; then
                log_info "${protein}/${method}: ${found} prediction files found"
            else
                log_error "${protein}/${method}: no prediction_*.csv found"
                ok=0
            fi
        done
    done

    if [ "$ok" -eq 0 ]; then
        log_error "Stage 2 verification failed"
        [ -n "$STATE_FILE" ] && update_stage_status "2" "failed" 1
        exit 1
    fi

    # Update state to completed
    [ -n "$STATE_FILE" ] && update_stage_status "2" "completed" 0

    log_stage_complete "Stage 2: Initial Screening"
}

main "$@"
