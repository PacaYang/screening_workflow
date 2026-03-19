#!/bin/bash
#
# Stage 5: Fine Screening
# Description: Run fine screening methods with dependencies (DiffDock depends on Vina, MD+PBSA depends on DiffDock)
# Sub-stages: 5a (submit AF3/Boltz2/Vina/RF), 5b (wait), 5c (submit DiffDock), 5d (wait), 5e (submit MD+PBSA), 5f (wait)
# Inputs: ${TASK_ROOT}/${PROTEIN}/fine_screening/*/input/
# Outputs: ${TASK_ROOT}/${PROTEIN}/fine_screening/*/output/
# Dependencies: Stage 4
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
SKIP_VINA=0
SKIP_DIFFDOCK=0
SKIP_MD_PBSA=0

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
        --skip-vina) SKIP_VINA=1; shift ;;
        --skip-diffdock) SKIP_DIFFDOCK=1; shift ;;
        --skip-md-pbsa) SKIP_MD_PBSA=1; shift ;;
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
    log_stage_start "Stage 5: Fine Screening"

    # Update state to running
    [ -n "$STATE_FILE" ] && update_stage_status "5" "running"

    # Export environment variables
    export_pipeline_env "$TASK_ROOT" "$PROTEINS"

    # 5a — Submit AF3, Boltz2, Vina, RoseTTAFold
    log_info "Stage 5a: Submitting AF3, Boltz2, Vina, RoseTTAFold jobs"
    [ "$SKIP_AF3" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../structure_prediction/run_af3_batch.sh"
    [ "$SKIP_BOLTZ2" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../structure_prediction/run_boltz2_batch.sh"
    [ "$SKIP_VINA" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../docking/run_vina_batch.sh"
    [ "$SKIP_ROSETTAFOLD" -eq 0 ] && run_or_dry "$DRY_RUN" bash "${SCRIPT_DIR}/../../structure_prediction/run_rosettafold_batch.sh"

    if [ "$DRY_RUN" -eq 1 ]; then
        log_info "[DRY RUN] Would wait for jobs and submit dependent stages"
        return 0
    fi

    # 5b — Wait for AF3, Boltz2, Vina, RoseTTAFold
    log_info "Stage 5b: Waiting for AF3, Boltz2, Vina, RoseTTAFold jobs"
    local job_files=()

    for protein in $PROTEINS; do
        [ "$SKIP_AF3" -eq 0 ] && [ -f "${TASK_ROOT}/${protein}/fine_screening/AF3/output/job_ids.txt" ] && \
            job_files+=("${TASK_ROOT}/${protein}/fine_screening/AF3/output/job_ids.txt")
        [ "$SKIP_BOLTZ2" -eq 0 ] && [ -f "${TASK_ROOT}/${protein}/fine_screening/Boltz2/output/job_ids.txt" ] && \
            job_files+=("${TASK_ROOT}/${protein}/fine_screening/Boltz2/output/job_ids.txt")
        [ "$SKIP_VINA" -eq 0 ] && [ -f "${TASK_ROOT}/${protein}/fine_screening/Vina/output/job_ids.txt" ] && \
            job_files+=("${TASK_ROOT}/${protein}/fine_screening/Vina/output/job_ids.txt")
        [ "$SKIP_ROSETTAFOLD" -eq 0 ] && [ -f "${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold/protein_ligand/job_ids.txt" ] && \
            job_files+=("${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold/protein_ligand/job_ids.txt")
    done

    if [ ${#job_files[@]} -gt 0 ]; then
        wait_for_jobs_with_progress "$POLL_INTERVAL" "${job_files[@]}"
    fi

    # 5c — Submit DiffDock (depends on Vina)
    if [ "$SKIP_DIFFDOCK" -eq 0 ]; then
        log_info "Stage 5c: Submitting DiffDock jobs"
        bash "${SCRIPT_DIR}/../../docking/run_diffdock_batch.sh"
    fi

    # 5d — Wait for DiffDock
    if [ "$SKIP_DIFFDOCK" -eq 0 ]; then
        log_info "Stage 5d: Waiting for DiffDock jobs"
        local dd_job_files=()
        for protein in $PROTEINS; do
            local jf="${TASK_ROOT}/${protein}/fine_screening/PBSA/DiffDock/output/job_ids.txt"
            [ -f "$jf" ] && dd_job_files+=("$jf")
        done
        if [ ${#dd_job_files[@]} -gt 0 ]; then
            wait_for_jobs_with_progress "$POLL_INTERVAL" "${dd_job_files[@]}"
        fi
    fi

    # 5e — Submit MD+PBSA (depends on DiffDock)
    if [ "$SKIP_MD_PBSA" -eq 0 ]; then
        log_info "Stage 5e: Submitting MD+PBSA jobs"
        bash "${SCRIPT_DIR}/../../md_pbsa/run_md_pbsa_batch.sh"
    fi

    # 5f — Wait for MD+PBSA
    if [ "$SKIP_MD_PBSA" -eq 0 ]; then
        log_info "Stage 5f: Waiting for MD+PBSA jobs"
        local pbsa_job_files=()
        for protein in $PROTEINS; do
            local jf="${TASK_ROOT}/${protein}/fine_screening/PBSA/PBSA/job_ids.txt"
            [ -f "$jf" ] && pbsa_job_files+=("$jf")
        done
        if [ ${#pbsa_job_files[@]} -gt 0 ]; then
            wait_for_jobs_with_progress "$POLL_INTERVAL" "${pbsa_job_files[@]}"
        fi
    fi

    # Update state to completed
    [ -n "$STATE_FILE" ] && update_stage_status "5" "completed" 0

    log_stage_complete "Stage 5: Fine Screening"
}

main "$@"
