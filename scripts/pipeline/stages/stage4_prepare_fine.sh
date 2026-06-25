#!/bin/bash
#
# Stage 4: Prepare Fine Screening
# Description: Write prefold inputs, run prefold jobs, write fine screening inputs
# Sub-stages: 4a (write prefold), 4b (submit prefold), 4c (wait), 4d (write fine)
# Inputs: ${TASK_ROOT}/${PROTEIN}/initial_screening/selected.csv
# Outputs: ${TASK_ROOT}/${PROTEIN}/fine_screening/*/input/
# Dependencies: Stage 3
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
SKIP_AF3=0
SKIP_BOLTZ2=0
SKIP_ROSETTAFOLD=0

# Template mode (single-stage, no prefold) and multi-ligand
TEMPLATE_MODE=0
N_LIGANDS=1

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --task-root) TASK_ROOT="$2"; shift 2 ;;
        --proteins) PROTEINS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --poll-interval) POLL_INTERVAL="$2"; shift 2 ;;
        --skip-af3) SKIP_AF3=1; shift ;;
        --skip-boltz2) SKIP_BOLTZ2=1; shift ;;
        --skip-rosettafold) SKIP_ROSETTAFOLD=1; shift ;;
        --template-mode) TEMPLATE_MODE=1; shift ;;
        --n-ligands) N_LIGANDS="$2"; shift 2 ;;
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
    log_stage_start "Stage 4: Prepare Fine Screening"

    # Update state to running
    [ -n "$STATE_FILE" ] && update_stage_status "4" "running"

    # Export environment variables
    export_pipeline_env "$TASK_ROOT" "$PROTEINS"
    # Make template settings visible to the input writers / batch scripts.
    export MASTER_TEMPLATE_MODE="$TEMPLATE_MODE"
    export MASTER_N_LIGANDS="$N_LIGANDS"

    # Template mode: single stage, no prefold/MSA. Write template inputs directly
    # from sequences.csv, then prepare the docking CSV splits as usual.
    if [ "$TEMPLATE_MODE" -eq 1 ]; then
        log_info "Stage 4 (template mode): writing template-based fine screening inputs"
        [ "$SKIP_AF3" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../input/run_write_af3_template_input.sh"
        [ "$SKIP_BOLTZ2" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../input/run_write_boltz2_template_input.sh"
        # RoseTTAFold template inputs are generated inside run_rosettafold_batch.sh (stage 5).

        if [ "$DRY_RUN" -eq 1 ]; then
            log_info "[DRY RUN] Would write template inputs and split docking CSV"
            return 0
        fi

        bash "${SCRIPT_DIR}/../../input/run_split_csv.sh"

        [ -n "$STATE_FILE" ] && update_stage_status "4" "completed" 0
        log_stage_complete "Stage 4: Prepare Fine Screening (template mode)"
        return 0
    fi

    # 4a — Write prefold inputs
    log_info "Stage 4a: Writing prefold inputs"
    [ "$SKIP_AF3" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../input/run_write_prefold_af3.sh"
    [ "$SKIP_BOLTZ2" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../input/run_write_prefold_boltz2.sh"

    # 4b — Submit prefold SLURM jobs
    log_info "Stage 4b: Submitting prefold jobs"
    [ "$SKIP_AF3" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../structure_prediction/prefold/run_prefold_af3.sh"
    [ "$SKIP_BOLTZ2" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../structure_prediction/prefold/run_prefold_boltz2.sh"
    [ "$SKIP_ROSETTAFOLD" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../structure_prediction/prefold/run_rosettafold_prefold.sh"

    if [ "$DRY_RUN" -eq 1 ]; then
        log_info "[DRY RUN] Would wait for prefold jobs and write fine screening inputs"
        return 0
    fi

    # 4c — Wait for prefold jobs
    log_info "Stage 4c: Waiting for prefold jobs"
    local job_files=()

    if [ "$SKIP_AF3" -eq 0 ]; then
        for protein in $PROTEINS; do
            local jf="${TASK_ROOT}/${protein}/fine_screening/AF3/prefold/job_ids.txt"
            [ -f "$jf" ] && job_files+=("$jf")
        done
    fi

    if [ "$SKIP_BOLTZ2" -eq 0 ]; then
        for protein in $PROTEINS; do
            local jf="${TASK_ROOT}/${protein}/fine_screening/Boltz2/prefold/job_ids.txt"
            [ -f "$jf" ] && job_files+=("$jf")
        done
    fi

    if [ "$SKIP_ROSETTAFOLD" -eq 0 ]; then
        for protein in $PROTEINS; do
            local jf="${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold/job_ids.txt"
            [ -f "$jf" ] && job_files+=("$jf")
        done
    fi

    if [ ${#job_files[@]} -gt 0 ]; then
        wait_for_jobs_with_progress "$POLL_INTERVAL" "${job_files[@]}"
    fi

    # 4d — Write fine screening inputs
    log_info "Stage 4d: Writing fine screening inputs"
    [ "$SKIP_AF3" -eq 0 ] && bash "${SCRIPT_DIR}/../../input/run_write_af3_input.sh"
    [ "$SKIP_BOLTZ2" -eq 0 ] && bash "${SCRIPT_DIR}/../../input/run_write_boltz2_input.sh"
    bash "${SCRIPT_DIR}/../../input/run_split_csv.sh"

    # Update state to completed
    [ -n "$STATE_FILE" ] && update_stage_status "4" "completed" 0

    log_stage_complete "Stage 4: Prepare Fine Screening"
}

main "$@"
