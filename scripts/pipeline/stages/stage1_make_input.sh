#!/bin/bash
#
# Stage 1: Make Input
# Description: Generate input files for initial screening
# Inputs: ${TASK_ROOT}/Input/ directory with compounds_smiles.csv, sequences.csv
# Outputs: ${TASK_ROOT}/${PROTEIN}/initial_screening/inputs/finish.token
# Dependencies: None
#

set -e

# Source libraries
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../lib/pipeline_utils.sh"
source "${SCRIPT_DIR}/../lib/state_manager.sh"
source "${SCRIPT_DIR}/../lib/logger.sh"

# Default values
TASK_ROOT=""
PROTEINS=""
DRY_RUN=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --task-root)
            TASK_ROOT="$2"
            shift 2
            ;;
        --proteins)
            PROTEINS="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --state-file)
            STATE_FILE="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$TASK_ROOT" ] || [ -z "$PROTEINS" ]; then
    log_error "Missing required arguments: --task-root and --proteins"
    exit 1
fi

# Main execution
main() {
    log_stage_start "Stage 1: Make Input"

    # Check if already complete
    if is_stage_complete "$TASK_ROOT" "$PROTEINS" "initial_screening/inputs/finish.token"; then
        log_info "Stage 1 already complete, skipping"
        [ -n "$STATE_FILE" ] && update_stage_status "1" "completed" 0
        return 0
    fi

    # Update state to running
    [ -n "$STATE_FILE" ] && update_stage_status "1" "running"

    # Export environment variables
    export_pipeline_env "$TASK_ROOT" "$PROTEINS"

    # Find and execute run_make_input.sh
    local make_input_script="${SCRIPT_DIR}/../../input/run_make_input.sh"
    check_script_exists "$make_input_script" "run_make_input.sh" || exit 1

    # Execute
    run_or_dry "$DRY_RUN" bash "$make_input_script"

    if [ "$DRY_RUN" -eq 1 ]; then
        log_info "[DRY RUN] Stage 1 would complete here"
        return 0
    fi

    # Validate outputs
    local ok=1
    for protein in $PROTEINS; do
        local token="${TASK_ROOT}/${protein}/initial_screening/inputs/finish.token"
        if [ -f "$token" ]; then
            log_info "${protein}: input ready"
        else
            log_error "${protein}: finish.token not found at ${token}"
            ok=0
        fi
    done

    if [ "$ok" -eq 0 ]; then
        log_error "Stage 1 verification failed"
        [ -n "$STATE_FILE" ] && update_stage_status "1" "failed" 1
        exit 1
    fi

    # Update state to completed
    [ -n "$STATE_FILE" ] && update_stage_status "1" "completed" 0

    log_stage_complete "Stage 1: Make Input"
}

main "$@"
