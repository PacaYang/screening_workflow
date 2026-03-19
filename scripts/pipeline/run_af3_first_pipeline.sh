#!/bin/bash
#
# AF3-First Pipeline (No Initial Screening)
# Flow:
#   1) Run all prefold jobs (AF3, Boltz2, RoseTTAFold)
#   2) Run AF3 on molecules from Input/compounds_smiles.csv
#   3) Select top fraction from AF3 results (pLDDT + PAE criteria)
#   4) Run Boltz2, RoseTTAFold, Vina on selected molecules
#   5) Collect final results
#

set -e

# Source libraries
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/logger.sh"
source "${SCRIPT_DIR}/lib/pipeline_utils.sh"
source "${SCRIPT_DIR}/lib/job_monitor.sh"

# ============================================================================
# Defaults
# ============================================================================

TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
PROTEINS="${MASTER_PROTEINS:-JAK1JH1}"
POLL_INTERVAL=300
DRY_RUN=0

AF3_PLDDT_THRESHOLD=70
TOP_FRACTION=0.3

MODEL_WEIGHTS_DIR="/home/yangl_pacagen_com/Applications/model_weights"

START_FROM=1
STOP_AFTER=6

COMMAND="run"

# Checkpoints
CHECKPOINT_DIR=""

# ============================================================================
# Helpers
# ============================================================================

show_help() {
    cat <<'HELPEOF'
AF3-First Pipeline (No Initial Screening)
=========================================

Usage: run_af3_first_pipeline.sh [OPTIONS] [COMMAND]

Commands:
  run              Run pipeline (default)
  status           Show AF3-first checkpoint status
  help             Show this help message

Options:
  --task-root DIR             Task root directory
  --proteins "P1 P2"          Space-separated protein list
  --poll-interval SEC         Job poll interval (default: 300)
  --af3-plddt-threshold N     AF3 selection pLDDT threshold (default: 70)
  --top-fraction F            Top fraction for AF3 selection (default: 0.3)
  --model-weights-dir DIR     Model weights root (default: /home/.../model_weights)
  --start-from STEP           Start from step 1-6 (default: 1)
  --stop-after STEP           Stop after step 1-6 (default: 6)
  --dry-run                   Show what would run

Notes:
  - Seed molecules are read from:
      ${TASK_ROOT}/Input/compounds_smiles.csv
  - This flow does not run initial screening methods.
HELPEOF
}

cleanup() {
    echo ""
    log_error "AF3-first pipeline interrupted"
    exit 130
}
trap cleanup SIGINT SIGTERM

parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --task-root) TASK_ROOT="$2"; shift 2 ;;
            --proteins) PROTEINS="$2"; shift 2 ;;
            --poll-interval) POLL_INTERVAL="$2"; shift 2 ;;
            --af3-plddt-threshold) AF3_PLDDT_THRESHOLD="$2"; shift 2 ;;
            --top-fraction) TOP_FRACTION="$2"; shift 2 ;;
            --model-weights-dir) MODEL_WEIGHTS_DIR="$2"; shift 2 ;;
            --start-from) START_FROM="$2"; shift 2 ;;
            --stop-after) STOP_AFTER="$2"; shift 2 ;;
            --dry-run) DRY_RUN=1; shift ;;
            run|status|help) COMMAND="$1"; shift ;;
            --help|-h) COMMAND="help"; shift ;;
            *) echo "Unknown option: $1"; show_help; exit 1 ;;
        esac
    done
}

validate_arguments() {
    if [ "$START_FROM" -lt 1 ] || [ "$START_FROM" -gt 6 ]; then
        log_error "--start-from must be 1-6"
        exit 1
    fi
    if [ "$STOP_AFTER" -lt 1 ] || [ "$STOP_AFTER" -gt 6 ]; then
        log_error "--stop-after must be 1-6"
        exit 1
    fi
    if [ "$START_FROM" -gt "$STOP_AFTER" ]; then
        log_error "--start-from cannot be greater than --stop-after"
        exit 1
    fi
}

setup_runtime() {
    CHECKPOINT_DIR="${TASK_ROOT}/.af3_first_state"
    mkdir -p "$CHECKPOINT_DIR"

    local log_dir="${TASK_ROOT}/.pipeline_logs"
    mkdir -p "$log_dir"
    LOG_FILE="${log_dir}/af3_first_$(date +%Y%m%d_%H%M%S).log"

    export_pipeline_env "$TASK_ROOT" "$PROTEINS"
    export MASTER_MODEL_WEIGHTS_DIR="$MODEL_WEIGHTS_DIR"
    export MASTER_AF3_WEIGHT_DIR="${MODEL_WEIGHTS_DIR}/AF3"
    export MASTER_AF3_DB_DIR="${MODEL_WEIGHTS_DIR}/af3_db"
}

step_name() {
    case "$1" in
        1) echo "Prefold All Methods" ;;
        2) echo "Bootstrap Seed Selection" ;;
        3) echo "AF3 Seed Run + Collect" ;;
        4) echo "AF3 Top-Fraction Select + Prepare Downstream" ;;
        5) echo "Run Downstream Screening" ;;
        6) echo "Collect Final Results" ;;
        *) echo "Unknown Step" ;;
    esac
}

step_token() {
    local step="$1"
    echo "${CHECKPOINT_DIR}/step${step}.done"
}

mark_step_done() {
    local step="$1"
    [ "$DRY_RUN" -eq 1 ] && return 0
    date -u +%Y-%m-%dT%H:%M:%SZ > "$(step_token "$step")"
}

is_step_done() {
    local step="$1"
    [ -f "$(step_token "$step")" ]
}

run_cmd() {
    if [ "$DRY_RUN" -eq 1 ]; then
        log_info "[DRY RUN] Would execute: $*"
        return 0
    fi
    "$@"
}

get_prefold_job_files() {
    local job_files=()
    local protein
    for protein in $PROTEINS; do
        [ -f "${TASK_ROOT}/${protein}/fine_screening/AF3/prefold/job_ids.txt" ] && \
            job_files+=("${TASK_ROOT}/${protein}/fine_screening/AF3/prefold/job_ids.txt")
        [ -f "${TASK_ROOT}/${protein}/fine_screening/Boltz2/prefold/job_ids.txt" ] && \
            job_files+=("${TASK_ROOT}/${protein}/fine_screening/Boltz2/prefold/job_ids.txt")
        [ -f "${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold/job_ids.txt" ] && \
            job_files+=("${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold/job_ids.txt")
    done
    echo "${job_files[@]}"
}

get_af3_job_files() {
    local job_files=()
    local protein
    for protein in $PROTEINS; do
        [ -f "${TASK_ROOT}/${protein}/fine_screening/AF3/output/job_ids.txt" ] && \
            job_files+=("${TASK_ROOT}/${protein}/fine_screening/AF3/output/job_ids.txt")
    done
    echo "${job_files[@]}"
}

get_downstream_job_files() {
    local job_files=()
    local protein
    for protein in $PROTEINS; do
        [ -f "${TASK_ROOT}/${protein}/fine_screening/Boltz2/output/job_ids.txt" ] && \
            job_files+=("${TASK_ROOT}/${protein}/fine_screening/Boltz2/output/job_ids.txt")
        [ -f "${TASK_ROOT}/${protein}/fine_screening/Vina/output/job_ids.txt" ] && \
            job_files+=("${TASK_ROOT}/${protein}/fine_screening/Vina/output/job_ids.txt")
        [ -f "${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold/protein_ligand/job_ids.txt" ] && \
            job_files+=("${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold/protein_ligand/job_ids.txt")
    done
    echo "${job_files[@]}"
}

run_step_1() {
    run_cmd bash "${SCRIPT_DIR}/../input/run_write_prefold_af3.sh"
    run_cmd bash "${SCRIPT_DIR}/../input/run_write_prefold_boltz2.sh"

    run_cmd bash "${SCRIPT_DIR}/../structure_prediction/prefold/run_prefold_af3.sh"
    run_cmd bash "${SCRIPT_DIR}/../structure_prediction/prefold/run_prefold_boltz2.sh"
    run_cmd bash "${SCRIPT_DIR}/../structure_prediction/prefold/run_rosettafold_prefold.sh"

    if [ "$DRY_RUN" -eq 1 ]; then
        return 0
    fi

    local job_files
    read -r -a job_files <<< "$(get_prefold_job_files)"
    if [ "${#job_files[@]}" -gt 0 ]; then
        wait_for_jobs_with_progress "$POLL_INTERVAL" "${job_files[@]}"
    else
        log_warn "No prefold job files found; continuing"
    fi

    local protein
    for protein in $PROTEINS; do
        local af3_token="${TASK_ROOT}/${protein}/fine_screening/AF3/prefold/prefold.done"
        local boltz_conf="${TASK_ROOT}/${protein}/fine_screening/Boltz2/prefold/boltz_results_${protein}/predictions/${protein}/confidence_${protein}_model_0.json"
        local rf_token="${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold/protein_folding/output/protein_fold.done"

        [ -f "$af3_token" ] || { log_error "Missing AF3 prefold token: ${af3_token}"; return 1; }
        [ -f "$boltz_conf" ] || { log_error "Missing Boltz2 prefold output: ${boltz_conf}"; return 1; }
        [ -f "$rf_token" ] || { log_error "Missing RoseTTAFold prefold token: ${rf_token}"; return 1; }
    done
}

run_step_2() {
    local seed_csv="${TASK_ROOT}/Input/compounds_smiles.csv"
    if [ "$DRY_RUN" -eq 0 ] && [ ! -f "$seed_csv" ]; then
        log_error "Seed CSV not found: ${seed_csv}"
        return 1
    fi

    local protein
    for protein in $PROTEINS; do
        local init_dir="${TASK_ROOT}/${protein}/initial_screening"
        run_cmd mkdir -p "$init_dir"
        run_cmd python "${SCRIPT_DIR}/../scoring/bootstrap_selected_from_seed.py" \
            --seed-csv "$seed_csv" \
            --output-selected "${init_dir}/selected.csv" \
            --output-seed-copy "${init_dir}/selected_seed.csv"

        if [ "$DRY_RUN" -eq 0 ]; then
            local n
            n=$(tail -n +2 "${init_dir}/selected.csv" | wc -l)
            log_info "${protein}: seeded ${n} compounds for AF3"
        fi
    done
}

run_step_3() {
    run_cmd bash "${SCRIPT_DIR}/../input/run_write_af3_input.sh"
    run_cmd bash "${SCRIPT_DIR}/../structure_prediction/run_af3_batch.sh"

    if [ "$DRY_RUN" -eq 1 ]; then
        return 0
    fi

    local job_files
    read -r -a job_files <<< "$(get_af3_job_files)"
    if [ "${#job_files[@]}" -gt 0 ]; then
        wait_for_jobs_with_progress "$POLL_INTERVAL" "${job_files[@]}"
    else
        log_warn "No AF3 job files found; continuing"
    fi

    local protein
    for protein in $PROTEINS; do
        local base="${TASK_ROOT}/${protein}/fine_screening/AF3"
        run_cmd python "${SCRIPT_DIR}/../scoring/af3_scores.py" \
            --af3-results-folder "${base}/output" \
            --output-dir "${base}"
        [ -f "${base}/summary.csv" ] || { log_error "${protein}: AF3 summary.csv missing"; return 1; }
    done
}

run_step_4() {
    local protein
    for protein in $PROTEINS; do
        local base="${TASK_ROOT}/${protein}"
        run_cmd python "${SCRIPT_DIR}/../scoring/select_af3_top_fraction.py" \
            --af3-summary "${base}/fine_screening/AF3/summary.csv" \
            --selected-output "${base}/initial_screening/selected.csv" \
            --top-output "${base}/fine_screening/AF3/selected_top_fraction.csv" \
            --plddt-threshold "${AF3_PLDDT_THRESHOLD}" \
            --top-fraction "${TOP_FRACTION}"

        if [ "$DRY_RUN" -eq 0 ]; then
            local n
            n=$(tail -n +2 "${base}/initial_screening/selected.csv" | wc -l)
            log_info "${protein}: selected ${n} compounds after AF3 filtering"
        fi
    done

    run_cmd bash "${SCRIPT_DIR}/../input/run_write_boltz2_input.sh"
    run_cmd bash "${SCRIPT_DIR}/../input/run_split_csv.sh"
}

run_step_5() {
    run_cmd bash "${SCRIPT_DIR}/../structure_prediction/run_boltz2_batch.sh"
    run_cmd bash "${SCRIPT_DIR}/../structure_prediction/run_rosettafold_batch.sh"
    run_cmd bash "${SCRIPT_DIR}/../docking/run_vina_batch.sh"

    if [ "$DRY_RUN" -eq 1 ]; then
        return 0
    fi

    local job_files
    read -r -a job_files <<< "$(get_downstream_job_files)"
    if [ "${#job_files[@]}" -gt 0 ]; then
        wait_for_jobs_with_progress "$POLL_INTERVAL" "${job_files[@]}"
    else
        log_warn "No downstream job files found; continuing"
    fi
}

run_step_6() {
    local protein
    for protein in $PROTEINS; do
        local base="${TASK_ROOT}/${protein}"

        if [ "$DRY_RUN" -eq 0 ] && [ ! -f "${base}/fine_screening/AF3/summary.csv" ]; then
            run_cmd python "${SCRIPT_DIR}/../scoring/af3_scores.py" \
                --af3-results-folder "${base}/fine_screening/AF3/output" \
                --output-dir "${base}/fine_screening/AF3"
        fi

        run_cmd python "${SCRIPT_DIR}/../scoring/boltz2_scores.py" \
            --boltz-results-folder "${base}/fine_screening/Boltz2/output" \
            --output-dir "${base}/fine_screening/Boltz2"

        run_cmd python "${SCRIPT_DIR}/../scoring/vina_scores.py" \
            --vina-results-folder "${base}/fine_screening/Vina/output" \
            --input-dir "${base}/fine_screening/Vina/input" \
            --output-dir "${base}/fine_screening/Vina"

        run_cmd conda run -n RFAA python "${SCRIPT_DIR}/../scoring/rosettafold_scores.py" \
            --rfaa-results-folder "${base}/fine_screening/RoseTTAFold/protein_ligand/output" \
            --protein-name "${protein}" \
            --output-dir "${base}/fine_screening/RoseTTAFold"

        if [ "$DRY_RUN" -eq 0 ]; then
            [ -f "${base}/fine_screening/AF3/summary.csv" ] || { log_error "${protein}: missing AF3 summary.csv"; return 1; }
            [ -f "${base}/fine_screening/Boltz2/summary.csv" ] || { log_error "${protein}: missing Boltz2 summary.csv"; return 1; }
            [ -f "${base}/fine_screening/Vina/results.csv" ] || { log_error "${protein}: missing Vina results.csv"; return 1; }
            [ -f "${base}/fine_screening/RoseTTAFold/summary.csv" ] || { log_error "${protein}: missing RoseTTAFold summary.csv"; return 1; }
        fi
    done
}

execute_step() {
    local step="$1"
    local step_label
    step_label="$(step_name "$step")"

    if [ "$step" -lt "$START_FROM" ] || [ "$step" -gt "$STOP_AFTER" ]; then
        log_info "Skipping step ${step} (${step_label}) due to range filter"
        return 0
    fi

    if is_step_done "$step"; then
        log_info "Skipping step ${step} (${step_label}) - checkpoint exists"
        return 0
    fi

    log_stage_start "Step ${step}: ${step_label}"
    "run_step_${step}"
    mark_step_done "$step"
    log_stage_complete "Step ${step}: ${step_label}"
}

run_pipeline() {
    setup_runtime

    log_info "Starting AF3-first pipeline"
    log_info "Task Root: ${TASK_ROOT}"
    log_info "Proteins: ${PROTEINS}"
    log_info "Steps: ${START_FROM} to ${STOP_AFTER}"
    log_info "Seed CSV: ${TASK_ROOT}/Input/compounds_smiles.csv"
    log_info "AF3 selection: pLDDT >= ${AF3_PLDDT_THRESHOLD}, PAE != -1, top fraction ${TOP_FRACTION}"
    log_info "DiffDock/MD-PBSA: excluded in this flow"

    local step
    for step in 1 2 3 4 5 6; do
        execute_step "$step"
    done

    log_info "AF3-first pipeline completed"
}

show_status() {
    CHECKPOINT_DIR="${TASK_ROOT}/.af3_first_state"
    echo "=========================================="
    echo "  AF3-First Pipeline Status"
    echo "  $(date)"
    echo "=========================================="
    echo "Task Root: ${TASK_ROOT}"
    echo "Proteins:  ${PROTEINS}"
    echo ""

    local step
    for step in 1 2 3 4 5 6; do
        local token
        token="$(step_token "$step")"
        if [ -f "$token" ]; then
            local ts
            ts=$(cat "$token")
            echo "Step ${step} ($(step_name "$step")): completed (${ts})"
        else
            echo "Step ${step} ($(step_name "$step")): pending"
        fi
    done
}

main() {
    parse_arguments "$@"
    validate_arguments

    case "$COMMAND" in
        run) run_pipeline ;;
        status) show_status ;;
        help) show_help ;;
        *)
            log_error "Unknown command: ${COMMAND}"
            show_help
            exit 1
            ;;
    esac
}

main "$@"
