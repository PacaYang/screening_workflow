#!/bin/bash
#
# Stage 6: Collect Results
# Description: Collect and compile final results from all fine screening methods
# Inputs: ${TASK_ROOT}/${PROTEIN}/fine_screening/*/output/
# Outputs: ${TASK_ROOT}/${PROTEIN}/fine_screening/*/summary.csv
# Dependencies: Stage 5
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
STREAMING_MODE=1  # Default to streaming
COLLECTION_INTERVAL=3600
MAX_ITERATIONS=0

# Skip flags
SKIP_AF3=0
SKIP_BOLTZ2=0
SKIP_VINA=0
SKIP_MD_PBSA=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --task-root) TASK_ROOT="$2"; shift 2 ;;
        --proteins) PROTEINS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --streaming-collection) STREAMING_MODE=1; shift ;;
        --batch-collection) STREAMING_MODE=0; shift ;;
        --collection-interval) COLLECTION_INTERVAL="$2"; shift 2 ;;
        --max-iterations) MAX_ITERATIONS="$2"; shift 2 ;;
        --skip-af3) SKIP_AF3=1; shift ;;
        --skip-boltz2) SKIP_BOLTZ2=1; shift ;;
        --skip-vina) SKIP_VINA=1; shift ;;
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

# Delegate to streaming or batch mode
if [ "$STREAMING_MODE" -eq 1 ]; then
    log_info "Using streaming collection mode"
    exec bash "${SCRIPT_DIR}/stage6_collect_results_streaming.sh" \
        --task-root "$TASK_ROOT" \
        --proteins "$PROTEINS" \
        $([ "$DRY_RUN" -eq 1 ] && echo "--dry-run") \
        --collection-interval "$COLLECTION_INTERVAL" \
        --max-iterations "$MAX_ITERATIONS" \
        $([ "$SKIP_AF3" -eq 1 ] && echo "--skip-af3") \
        $([ "$SKIP_BOLTZ2" -eq 1 ] && echo "--skip-boltz2") \
        $([ "$SKIP_VINA" -eq 1 ] && echo "--skip-vina") \
        $([ "$SKIP_MD_PBSA" -eq 1 ] && echo "--skip-md-pbsa") \
        $([ -n "$STATE_FILE" ] && echo "--state-file $STATE_FILE")
fi

# Batch mode (original behavior)
log_info "Using batch collection mode"

# Main execution
main() {
    log_stage_start "Stage 6: Collect Results"

    # Update state to running
    [ -n "$STATE_FILE" ] && update_stage_status "6" "running"

    if [ "$DRY_RUN" -eq 1 ]; then
        log_info "[DRY RUN] Would collect scores for AF3, Boltz2, Vina, PBSA"
        return 0
    fi

    for protein in $PROTEINS; do
        log_info "Collecting results for ${protein}"
        local base="${TASK_ROOT}/${protein}"
        local selected="${base}/initial_screening/selected.csv"

        # AF3 scores
        if [ "$SKIP_AF3" -eq 0 ]; then
            local af3_summary="${base}/fine_screening/AF3/summary.csv"
            if [ -f "$af3_summary" ]; then
                log_info "${protein}/AF3: summary.csv already exists, skipping"
            else
                log_info "${protein}/AF3: collecting scores"
                local af3_out="${base}/fine_screening/AF3/output"
                python "${SCRIPT_DIR}/../../scoring/af3_scores.py" \
                    --af3-results-folder "${af3_out}" \
                    --output-dir "${base}/fine_screening/AF3"
            fi
        fi

        # Boltz2 scores
        if [ "$SKIP_BOLTZ2" -eq 0 ]; then
            local boltz2_summary="${base}/fine_screening/Boltz2/summary.csv"
            if [ -f "$boltz2_summary" ]; then
                log_info "${protein}/Boltz2: summary.csv already exists, skipping"
            else
                log_info "${protein}/Boltz2: collecting scores"
                local boltz_out="${base}/fine_screening/Boltz2/output"
                python "${SCRIPT_DIR}/../../scoring/boltz2_scores.py" \
                    --boltz-results-folder "${boltz_out}" \
                    --output-dir "${base}/fine_screening/Boltz2"
            fi
        fi

        # Vina scores
        if [ "$SKIP_VINA" -eq 0 ]; then
            local vina_results="${base}/fine_screening/Vina/results.csv"
            if [ -f "$vina_results" ]; then
                log_info "${protein}/Vina: results.csv already exists, skipping"
            else
                log_info "${protein}/Vina: collecting scores"
                python "${SCRIPT_DIR}/../../scoring/vina_scores.py" \
                    --vina-results-folder "${base}/fine_screening/Vina/output" \
                    --output-dir "${base}/fine_screening/Vina" \
                    --input-dir "${base}/fine_screening/Vina/input"
            fi
        fi

        # PBSA scores
        if [ "$SKIP_MD_PBSA" -eq 0 ]; then
            local pbsa_summary="${base}/fine_screening/PBSA/summary.csv"
            if [ -f "$pbsa_summary" ]; then
                log_info "${protein}/PBSA: summary.csv already exists, skipping"
            else
                local pbsa_dir="${base}/fine_screening/PBSA/PBSA/PBSA"
                local pbsa_outdir="${base}/fine_screening/PBSA"
                log_info "${protein}/PBSA: extracting results"
                bash "${SCRIPT_DIR}/../../md_pbsa/pbsa/pbsa_extract_results.sh" "$pbsa_dir" "$pbsa_outdir"
                log_info "${protein}/PBSA: mapping SMILES"
                python "${SCRIPT_DIR}/../../md_pbsa/pbsa/mapping_smiles.py" \
                    --collected "${pbsa_outdir}/tmp.csv" \
                    --smiles_csv "$selected" \
                    --outdir "$pbsa_outdir"
            fi
        fi
    done

    # Update state to completed
    [ -n "$STATE_FILE" ] && update_stage_status "6" "completed" 0

    log_stage_complete "Stage 6: Collect Results"
}

main "$@"