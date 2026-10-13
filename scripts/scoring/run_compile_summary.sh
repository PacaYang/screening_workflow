#!/bin/bash
#
# Standalone script for Snakemake checkpoint: compile_per_protein_summary
# Aggregates initial screening results and selects top compounds.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
EXE="/home/yangl_pacagen_com/screening_workflow/scripts/scoring/init_select_top.py"
TARGET_N="${MASTER_TARGET_N:-3000}"
# Priority selection: fill a priority tier with an even DrugLAMP+HMSA split
# before topping up from the other methods. Enabled when MASTER_PRIORITY_MODE=1.
PRIORITY_MODE="${MASTER_PRIORITY_MODE:-0}"
PRIORITY_N="${MASTER_PRIORITY_N:-4000}"

# ============================================================================
# Functions
# ============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

get_proteins() {
    if [ -n "$MASTER_PROTEINS" ]; then
        echo "$MASTER_PROTEINS"
    else
        echo "JAK1JH1"
    fi
}

# ============================================================================
# Main
# ============================================================================

log_info "Starting compile_per_protein_summary"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    INIT_DIR="${TASK_ROOT}/${PROTEIN}/initial_screening"
    SUMMARY="${INIT_DIR}/summary.csv"
    SELECTED="${INIT_DIR}/selected.csv"

    if [ -f "$SELECTED" ]; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    # Check prerequisites. Only require methods that were actually run; skipped
    # methods (via MASTER_SKIP_<METHOD>=1) are excluded so their absence does not
    # block selection. init_select_top.py tolerates missing optional dirs.
    REQUIRED_METHODS=(GraphDTA HMSA ColdDTA)
    [ "${MASTER_SKIP_DRUGLAMP:-0}" -ne 1 ] && REQUIRED_METHODS+=(DrugLAMP)
    [ "${MASTER_SKIP_CONPLEX:-0}" -ne 1 ] && REQUIRED_METHODS+=(ConPLex)
    [ "${MASTER_SKIP_GRAPHDTA:-0}" -eq 1 ] && REQUIRED_METHODS=("${REQUIRED_METHODS[@]/GraphDTA}")
    [ "${MASTER_SKIP_HMSA:-0}" -eq 1 ] && REQUIRED_METHODS=("${REQUIRED_METHODS[@]/HMSA}")
    [ "${MASTER_SKIP_COLDDTA:-0}" -eq 1 ] && REQUIRED_METHODS=("${REQUIRED_METHODS[@]/ColdDTA}")

    MISSING=0
    for METHOD_DIR in "${REQUIRED_METHODS[@]}"; do
        [ -z "$METHOD_DIR" ] && continue
        if ! ls "${INIT_DIR}/${METHOD_DIR}"/prediction_*.csv &>/dev/null; then
            log_error "No prediction CSVs found in ${INIT_DIR}/${METHOD_DIR}/"
            MISSING=1
        fi
    done
    if [ "$MISSING" -eq 1 ]; then
        log_error "Skipping ${PROTEIN} due to missing prerequisites"
        continue
    fi

    conda activate HMSA_test

    mkdir -p "$INIT_DIR"
    PRIORITY_ARGS=()
    if [ "$PRIORITY_MODE" -eq 1 ]; then
        PRIORITY_ARGS=(--priority-mode --priority-n "$PRIORITY_N")
    fi
    python "$EXE" \
        --graphdta-dir "${INIT_DIR}/GraphDTA" \
        --hmsa-dir "${INIT_DIR}/HMSA" \
        --colddta-dir "${INIT_DIR}/ColdDTA" \
        --druglamp-dir "${INIT_DIR}/DrugLAMP" \
        --conplex-dir "${INIT_DIR}/ConPLex" \
        --target-n $TARGET_N \
        "${PRIORITY_ARGS[@]}" \
        --summary "$SUMMARY" \
        --selected "$SELECTED"

    log_info "Completed compile_per_protein_summary for ${PROTEIN}"
done

log_info "Done"
