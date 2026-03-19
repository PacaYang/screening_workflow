#!/bin/bash
#
# Master Controller for Modular Pipeline
# Orchestrates all 6 stages with dependency management, state tracking, and parallel execution
#

set -e

# Source libraries
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/logger.sh"
source "${SCRIPT_DIR}/lib/pipeline_utils.sh"
source "${SCRIPT_DIR}/lib/state_manager.sh"
source "${SCRIPT_DIR}/lib/job_monitor.sh"
source "${SCRIPT_DIR}/lib/config_loader.sh"

# ============================================================================
# Configuration Defaults
# ============================================================================

TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
PROTEINS="${MASTER_PROTEINS:-JAK1JH1}"

# Stage control
START_FROM=1
STOP_AFTER=6
POLL_INTERVAL=300
DRY_RUN=0

# Skip flags — initial screening
SKIP_GRAPHDTA=0
SKIP_HMSA=0
SKIP_COLDDTA=0
SKIP_DRUGLAMP=0
SKIP_CONPLEX=0

# Skip flags — fine screening
SKIP_AF3=0
SKIP_BOLTZ2=0
SKIP_ROSETTAFOLD=0
SKIP_VINA=0
SKIP_DIFFDOCK=0
SKIP_MD_PBSA=0

# Configuration
TARGET_N=3000
CONFIG_FILE=""
RESUME_ID=""
COMMAND="run"

# Streaming collection configuration
STREAMING_COLLECTION=1  # Default to streaming
COLLECTION_INTERVAL=3600
MAX_COLLECTION_ITERATIONS=0

# Model weights configuration
MODEL_WEIGHTS_DIR="/home/yangl_pacagen_com/Applications/model_weights"

# ============================================================================
# Trap handler for interruptions
# ============================================================================

cleanup() {
    echo ""
    log_error "Pipeline interrupted"
    if [ -n "$STATE_FILE" ]; then
        mark_pipeline_complete "interrupted"
        log_info "Pipeline state saved. Resume with: --resume-id ${PIPELINE_ID}"
    fi
    exit 130
}
trap cleanup SIGINT SIGTERM

# ============================================================================
# Help
# ============================================================================

show_help() {
    cat <<'HELPEOF'
Full Automation Pipeline for Screening Workflow (Modular Version)
==================================================================

Usage: run_full_pipeline.sh [OPTIONS] [COMMAND]

Commands:
  run              Run pipeline (default)
  status           Show pipeline status
  resume           Resume interrupted pipeline
  help             Show this help message

Options:
  --task-root DIR          Task root directory
  --proteins "P1 P2"       Space-separated protein list
  --config FILE            Load configuration from YAML file
  --start-from STAGE       Resume from stage 1-6 (default: 1)
  --stop-after STAGE       Stop after stage 1-6 (default: 6)
  --poll-interval SEC      Job poll interval (default: 300)
  --skip-<method>          Skip specific methods
  --target-n N             Compounds for fine screening (default: 3000)
  --streaming-collection   Use streaming results collection (default)
  --batch-collection       Use batch results collection
  --collection-interval N  Seconds between collection runs (default: 3600)
  --max-collection-iter N  Max collection iterations (default: 0=infinite)
  --dry-run                Show what would run
  --resume-id ID           Resume specific pipeline run
HELPEOF
}

# ============================================================================
# Status command
# ============================================================================

show_status() {
    if [ -n "$RESUME_ID" ]; then
        load_pipeline_state "$RESUME_ID" "$TASK_ROOT"
    else
        local latest_id
        latest_id=$(get_latest_pipeline_id "$TASK_ROOT")
        if [ -n "$latest_id" ]; then
            load_pipeline_state "$latest_id" "$TASK_ROOT"
        else
            log_error "No pipeline state found"
            return 1
        fi
    fi

    echo "=========================================="
    echo "  Pipeline Status Report"
    echo "  $(date)"
    echo "=========================================="
    echo ""
    echo "Pipeline ID: ${PIPELINE_ID}"
    echo "Task Root:   ${TASK_ROOT}"
    echo "Proteins:    $(jq -r '.proteins | join(" ")' "$STATE_FILE")"
    echo "Status:      $(jq -r '.status' "$STATE_FILE")"
    echo "Progress:    $(get_pipeline_progress)%"
    echo ""

    for stage in {1..6}; do
        local status
        status=$(get_stage_status "$stage")
        echo "Stage ${stage}: ${status}"
    done
}

# ============================================================================
# Argument Parsing
# ============================================================================

parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --task-root) TASK_ROOT="$2"; shift 2 ;;
            --proteins) PROTEINS="$2"; shift 2 ;;
            --config) CONFIG_FILE="$2"; shift 2 ;;
            --start-from) START_FROM="$2"; shift 2 ;;
            --stop-after) STOP_AFTER="$2"; shift 2 ;;
            --poll-interval) POLL_INTERVAL="$2"; shift 2 ;;
            --skip-graphdta) SKIP_GRAPHDTA=1; shift ;;
            --skip-hmsa) SKIP_HMSA=1; shift ;;
            --skip-colddta) SKIP_COLDDTA=1; shift ;;
            --skip-druglamp) SKIP_DRUGLAMP=1; shift ;;
            --skip-conplex) SKIP_CONPLEX=1; shift ;;
            --skip-af3) SKIP_AF3=1; shift ;;
            --skip-boltz2) SKIP_BOLTZ2=1; shift ;;
            --skip-rosettafold) SKIP_ROSETTAFOLD=1; shift ;;
            --skip-vina) SKIP_VINA=1; shift ;;
            --skip-diffdock) SKIP_DIFFDOCK=1; shift ;;
            --skip-md-pbsa) SKIP_MD_PBSA=1; shift ;;
            --target-n) TARGET_N="$2"; shift 2 ;;
            --streaming-collection) STREAMING_COLLECTION=1; shift ;;
            --batch-collection) STREAMING_COLLECTION=0; shift ;;
            --collection-interval) COLLECTION_INTERVAL="$2"; shift 2 ;;
            --max-collection-iter) MAX_COLLECTION_ITERATIONS="$2"; shift 2 ;;
            --dry-run) DRY_RUN=1; shift ;;
            --resume-id) RESUME_ID="$2"; shift 2 ;;
            --model-weights-dir) MODEL_WEIGHTS_DIR="$2"; shift 2 ;;
            run|status|resume|help) COMMAND="$1"; shift ;;
            all) COMMAND="run"; shift ;;
            --help|-h) COMMAND="help"; shift ;;
            *) echo "Unknown option: $1"; show_help; exit 1 ;;
        esac
    done
}

# ============================================================================
# Stage Execution
# ============================================================================

execute_stage() {
    local stage="$1"
    local stage_script="${SCRIPT_DIR}/stages/stage${stage}_*.sh"
    stage_script=$(ls $stage_script 2>/dev/null | head -1)

    if [ ! -f "$stage_script" ]; then
        log_error "Stage script not found for stage ${stage}"
        return 1
    fi

    log_info "Executing Stage ${stage}"

    # Build arguments
    local args=(
        --task-root "$TASK_ROOT"
        --proteins "$PROTEINS"
        --state-file "$STATE_FILE"
    )

    [ "$DRY_RUN" -eq 1 ] && args+=(--dry-run)

    # Add poll-interval for stages that need it (2, 4, 5)
    if [ "$stage" -eq 2 ] || [ "$stage" -eq 4 ] || [ "$stage" -eq 5 ]; then
        args+=(--poll-interval "$POLL_INTERVAL")
    fi

    # Add stage-specific skip flags
    case $stage in
        2)
            [ "$SKIP_GRAPHDTA" -eq 1 ] && args+=(--skip-graphdta)
            [ "$SKIP_HMSA" -eq 1 ] && args+=(--skip-hmsa)
            [ "$SKIP_COLDDTA" -eq 1 ] && args+=(--skip-colddta)
            [ "$SKIP_DRUGLAMP" -eq 1 ] && args+=(--skip-druglamp)
            [ "$SKIP_CONPLEX" -eq 1 ] && args+=(--skip-conplex)
            ;;
        4)
            [ "$SKIP_AF3" -eq 1 ] && args+=(--skip-af3)
            [ "$SKIP_BOLTZ2" -eq 1 ] && args+=(--skip-boltz2)
            [ "$SKIP_ROSETTAFOLD" -eq 1 ] && args+=(--skip-rosettafold)
            ;;
        5)
            [ "$SKIP_AF3" -eq 1 ] && args+=(--skip-af3)
            [ "$SKIP_BOLTZ2" -eq 1 ] && args+=(--skip-boltz2)
            [ "$SKIP_ROSETTAFOLD" -eq 1 ] && args+=(--skip-rosettafold)
            [ "$SKIP_VINA" -eq 1 ] && args+=(--skip-vina)
            [ "$SKIP_DIFFDOCK" -eq 1 ] && args+=(--skip-diffdock)
            [ "$SKIP_MD_PBSA" -eq 1 ] && args+=(--skip-md-pbsa)
            ;;
        6)
            [ "$SKIP_AF3" -eq 1 ] && args+=(--skip-af3)
            [ "$SKIP_BOLTZ2" -eq 1 ] && args+=(--skip-boltz2)
            [ "$SKIP_VINA" -eq 1 ] && args+=(--skip-vina)
            [ "$SKIP_ROSETTAFOLD" -eq 1 ] && args+=(--skip-rosettafold)
            [ "$SKIP_MD_PBSA" -eq 1 ] && args+=(--skip-md-pbsa)
            # Add streaming collection parameters
            [ "$STREAMING_COLLECTION" -eq 1 ] && args+=(--streaming-collection) || args+=(--batch-collection)
            args+=(--collection-interval "$COLLECTION_INTERVAL")
            args+=(--max-iterations "$MAX_COLLECTION_ITERATIONS")
            ;;
    esac

    # Execute stage script
    bash "$stage_script" "${args[@]}"
    return $?
}

# ============================================================================
# Main Pipeline Execution
# ============================================================================

run_pipeline() {
    log_info "Starting Pipeline"
    log_info "Task Root: ${TASK_ROOT}"
    log_info "Proteins: ${PROTEINS}"
    log_info "Stages: ${START_FROM} to ${STOP_AFTER}"

    # Initialize or resume pipeline state
    if [ -n "$RESUME_ID" ]; then
        log_info "Resuming pipeline: ${RESUME_ID}"
        load_pipeline_state "$RESUME_ID" "$TASK_ROOT"
    else
        init_pipeline_state "$TASK_ROOT" "$PROTEINS"
    fi

    # Setup logging
    local log_dir="${TASK_ROOT}/.pipeline_logs"
    mkdir -p "$log_dir"
    LOG_FILE="${log_dir}/${PIPELINE_ID}.log"
    log_info "Log file: ${LOG_FILE}"

    # Export environment variables
    export_pipeline_env "$TASK_ROOT" "$PROTEINS"
    export MASTER_TARGET_N="$TARGET_N"
    export MASTER_MODEL_WEIGHTS_DIR="$MODEL_WEIGHTS_DIR"
    export MASTER_AF3_WEIGHT_DIR="${MODEL_WEIGHTS_DIR}/AF3"
    export MASTER_AF3_DB_DIR="${MODEL_WEIGHTS_DIR}/af3_db"

    # Execute stages
    for stage in $(seq $START_FROM $STOP_AFTER); do
        log_info "========== Stage ${stage} =========="

        if execute_stage "$stage"; then
            log_info "Stage ${stage} completed successfully"
        else
            log_error "Stage ${stage} failed"
            mark_pipeline_complete "failed"
            exit 1
        fi
    done

    # Mark pipeline as complete
    mark_pipeline_complete "completed"
    log_info "Pipeline completed successfully!"
    log_info "Pipeline ID: ${PIPELINE_ID}"
}

# ============================================================================
# Main Entry Point
# ============================================================================

main() {
    parse_arguments "$@"

    case "$COMMAND" in
        run)
            run_pipeline
            ;;
        status)
            show_status
            ;;
        resume)
            if [ -z "$RESUME_ID" ]; then
                RESUME_ID=$(get_latest_pipeline_id "$TASK_ROOT")
                if [ -z "$RESUME_ID" ]; then
                    log_error "No pipeline to resume"
                    exit 1
                fi
            fi
            run_pipeline
            ;;
        help)
            show_help
            ;;
        *)
            log_error "Unknown command: ${COMMAND}"
            show_help
            exit 1
            ;;
    esac
}

main "$@"
