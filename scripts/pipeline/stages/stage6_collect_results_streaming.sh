#!/bin/bash
#
# Stage 6: Streaming Results Collection
# Description: Monitor and collect results incrementally as they complete
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
source "${SCRIPT_DIR}/../lib/collection_state.sh"

# Default values
TASK_ROOT=""
PROTEINS=""
DRY_RUN=0
COLLECTION_INTERVAL=3600  # 1 hour
MAX_ITERATIONS=0          # 0 = infinite
VALIDATE_RESULTS=1

# Skip flags
SKIP_AF3=0
SKIP_BOLTZ2=0
SKIP_VINA=0
SKIP_ROSETTAFOLD=0
SKIP_MD_PBSA=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --task-root) TASK_ROOT="$2"; shift 2 ;;
        --proteins) PROTEINS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --collection-interval) COLLECTION_INTERVAL="$2"; shift 2 ;;
        --max-iterations) MAX_ITERATIONS="$2"; shift 2 ;;
        --validate-results) VALIDATE_RESULTS=1; shift ;;
        --no-validate-results) VALIDATE_RESULTS=0; shift ;;
        --skip-af3) SKIP_AF3=1; shift ;;
        --skip-boltz2) SKIP_BOLTZ2=1; shift ;;
        --skip-vina) SKIP_VINA=1; shift ;;
        --skip-rosettafold) SKIP_ROSETTAFOLD=1; shift ;;
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

# Collect method incrementally
# Args: method, protein, base_dir, expected_count
# Returns: 0 if complete, 1 if still has uncollected
collect_method_incremental() {
    local method="$1"
    local protein="$2"
    local base="$3"
    local expected_count="$4"

    local state_dir="${base}/fine_screening/.collection_state"
    mkdir -p "$state_dir"
    local state_file="${state_dir}/${method}.json"

    # Initialize state if not exists
    if [ ! -f "$state_file" ]; then
        init_collection_state "$protein" "$method" "$expected_count" "$state_file"
    fi

    # Get collection statistics
    local stats
    stats=$(get_collection_stats "$state_file")
    local collected_count
    collected_count=$(echo "$stats" | jq -r '.collected')
    local uncollected=$((expected_count - collected_count))

    if [ $uncollected -eq 0 ]; then
        log_info "${protein}/${method}: All ${expected_count} compounds collected"
        return 0  # Complete
    fi

    log_info "${protein}/${method}: Collecting ${uncollected} uncollected compounds (${collected_count}/${expected_count} done)"

    # Call scoring script in incremental mode
    case "$method" in
        AF3)
            python "${SCRIPT_DIR}/../../scoring/af3_scores.py" \
                --af3-results-folder "${base}/fine_screening/AF3/output" \
                --output-dir "${base}/fine_screening/AF3" \
                --incremental \
                --state-file "$state_file" \
                --append
            ;;
        Boltz2)
            python "${SCRIPT_DIR}/../../scoring/boltz2_scores.py" \
                --boltz-results-folder "${base}/fine_screening/Boltz2/output" \
                --output-dir "${base}/fine_screening/Boltz2" \
                --incremental \
                --state-file "$state_file" \
                --append
            ;;
        Vina)
            python "${SCRIPT_DIR}/../../scoring/vina_scores.py" \
                --vina-results-folder "${base}/fine_screening/Vina/output" \
                --input-dir "${base}/fine_screening/Vina/input" \
                --output-dir "${base}/fine_screening/Vina"
            return 0
            ;;
        RoseTTAFold)
            conda run -n RFAA python "${SCRIPT_DIR}/../../scoring/rosettafold_scores.py" \
                --rfaa-results-folder "${base}/fine_screening/RoseTTAFold/protein_ligand/output" \
                --protein-name "$protein" \
                --output-dir "${base}/fine_screening/RoseTTAFold"
            return 0
            ;;
        PBSA)
            log_warn "${protein}/PBSA: Incremental collection not yet implemented, skipping"
            return 0
            ;;
    esac

    # Validate results
    if [ "$VALIDATE_RESULTS" -eq 1 ]; then
        validate_method_results "$method" "$protein" "$base" "$expected_count"
    fi

    # Re-check collected count after scoring run
    local new_stats
    new_stats=$(get_collection_stats "$state_file")
    local new_collected
    new_collected=$(echo "$new_stats" | jq -r '.collected')
    if [ "$new_collected" -ge "$expected_count" ]; then
        log_info "${protein}/${method}: All ${expected_count} compounds collected"
        return 0
    fi

    return 1  # Still has uncollected
}

# Validate method results
validate_method_results() {
    local method="$1"
    local protein="$2"
    local base="$3"
    local expected_count="$4"

    local summary_csv
    case "$method" in
        AF3) summary_csv="${base}/fine_screening/AF3/summary.csv" ;;
        Boltz2) summary_csv="${base}/fine_screening/Boltz2/summary.csv" ;;
        Vina) summary_csv="${base}/fine_screening/Vina/results.csv" ;;
        RoseTTAFold) summary_csv="${base}/fine_screening/RoseTTAFold/summary.csv" ;;
        PBSA) summary_csv="${base}/fine_screening/PBSA/summary.csv" ;;
    esac

    if [ ! -f "$summary_csv" ]; then
        log_warn "${protein}/${method}: Summary CSV not found"
        return
    fi

    # Check if empty
    local row_count
    row_count=$(tail -n +2 "$summary_csv" | wc -l)

    if [ $row_count -eq 0 ]; then
        log_error "${protein}/${method}: Summary CSV is EMPTY (0 rows)"
    elif [ $row_count -lt $((expected_count / 2)) ]; then
        local pct=$((row_count * 100 / expected_count))
        log_warn "${protein}/${method}: Only ${row_count}/${expected_count} compounds collected (${pct}%)"
    else
        log_info "${protein}/${method}: ${row_count}/${expected_count} compounds collected"
    fi
}

# Main execution
main() {
    log_stage_start "Stage 6: Streaming Results Collection"

    # Update state to running
    [ -n "$STATE_FILE" ] && update_stage_status "6" "running"

    if [ "$DRY_RUN" -eq 1 ]; then
        log_info "[DRY RUN] Would collect scores incrementally for AF3, Boltz2, Vina, PBSA"
        return 0
    fi

    local iteration=0
    local all_complete=0

    while [ $all_complete -eq 0 ]; do
        iteration=$((iteration + 1))
        log_info "Collection iteration ${iteration}"

        all_complete=1  # Assume complete unless we find uncollected

        for protein in $PROTEINS; do
            local base="${TASK_ROOT}/${protein}"
            local selected="${base}/initial_screening/selected.csv"

            # Get expected compound count
            if [ ! -f "$selected" ]; then
                log_warn "${protein}: selected.csv not found, skipping"
                continue
            fi

            local expected_count
            expected_count=$(tail -n +2 "$selected" | wc -l)

            # Collect each method incrementally
            if [ "$SKIP_AF3" -eq 0 ]; then
                if ! collect_method_incremental "AF3" "$protein" "$base" "$expected_count"; then
                    all_complete=0  # Still has uncollected
                fi
            fi

            if [ "$SKIP_BOLTZ2" -eq 0 ]; then
                if ! collect_method_incremental "Boltz2" "$protein" "$base" "$expected_count"; then
                    all_complete=0
                fi
            fi

            if [ "$SKIP_VINA" -eq 0 ]; then
                if ! collect_method_incremental "Vina" "$protein" "$base" "$expected_count"; then
                    all_complete=0
                fi
            fi

            if [ "$SKIP_ROSETTAFOLD" -eq 0 ]; then
                if ! collect_method_incremental "RoseTTAFold" "$protein" "$base" "$expected_count"; then
                    all_complete=0
                fi
            fi

            if [ "$SKIP_MD_PBSA" -eq 0 ]; then
                if ! collect_method_incremental "PBSA" "$protein" "$base" "$expected_count"; then
                    all_complete=0
                fi
            fi
        done

        # Check if we should continue
        if [ $all_complete -eq 1 ]; then
            log_info "All results collected"
            break
        fi

        if [ $MAX_ITERATIONS -gt 0 ] && [ $iteration -ge $MAX_ITERATIONS ]; then
            log_warn "Reached maximum iterations (${MAX_ITERATIONS}), stopping"
            break
        fi

        # Wait before next collection
        log_info "Waiting ${COLLECTION_INTERVAL}s before next collection..."
        sleep "$COLLECTION_INTERVAL"
    done

    # Update state to completed
    [ -n "$STATE_FILE" ] && update_stage_status "6" "completed" 0

    log_stage_complete "Stage 6: Streaming Results Collection"
}

main "$@"
