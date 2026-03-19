#!/bin/bash
#
# Stage 3: Compile Summary
# Description: Compile initial screening results and select top compounds
# Inputs: ${TASK_ROOT}/${PROTEIN}/initial_screening/${METHOD}/prediction_*.csv
# Outputs: ${TASK_ROOT}/${PROTEIN}/initial_screening/selected.csv
# Dependencies: Stage 2
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
        --task-root) TASK_ROOT="$2"; shift 2 ;;
        --proteins) PROTEINS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
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
    log_stage_start "Stage 3: Compile Summary"

    # Check if already complete
    if is_stage_complete "$TASK_ROOT" "$PROTEINS" "initial_screening/selected.csv"; then
        log_info "Stage 3 already complete, skipping"
        [ -n "$STATE_FILE" ] && update_stage_status "3" "completed" 0
        return 0
    fi

    # Update state to running
    [ -n "$STATE_FILE" ] && update_stage_status "3" "running"

    # Export environment variables
    export_pipeline_env "$TASK_ROOT" "$PROTEINS"

    # Find and execute run_compile_summary.sh
    local compile_script="${SCRIPT_DIR}/../../scoring/run_compile_summary.sh"
    check_script_exists "$compile_script" "run_compile_summary.sh" || exit 1

    # Execute
    run_or_dry "$DRY_RUN" bash "$compile_script"

    if [ "$DRY_RUN" -eq 1 ]; then
        log_info "[DRY RUN] Stage 3 would complete here"
        return 0
    fi

    # Validate outputs
    local ok=1
    for protein in $PROTEINS; do
        local sel="${TASK_ROOT}/${protein}/initial_screening/selected.csv"
        if [ -f "$sel" ]; then
            local n
            n=$(tail -n +2 "$sel" | wc -l)
            log_info "${protein}: ${n} compounds selected"
        else
            log_error "${protein}: selected.csv not found"
            ok=0
        fi
    done

    if [ "$ok" -eq 0 ]; then
        log_error "Stage 3 verification failed"
        [ -n "$STATE_FILE" ] && update_stage_status "3" "failed" 1
        exit 1
    fi

    # Update state to completed
    [ -n "$STATE_FILE" ] && update_stage_status "3" "completed" 0

    log_stage_complete "Stage 3: Compile Summary"
}

main "$@"