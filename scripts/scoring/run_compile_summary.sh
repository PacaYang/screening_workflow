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

    # Check prerequisites
    MISSING=0
    for METHOD_DIR in GraphDTA HMSA ColdDTA DrugLAMP ConPLex; do
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
    python "$EXE" \
        --graphdta-dir "${INIT_DIR}/GraphDTA" \
        --hmsa-dir "${INIT_DIR}/HMSA" \
        --colddta-dir "${INIT_DIR}/ColdDTA" \
        --druglamp-dir "${INIT_DIR}/DrugLAMP" \
        --conplex-dir "${INIT_DIR}/ConPLex" \
        --target-n $TARGET_N \
        --summary "$SUMMARY" \
        --selected "$SELECTED"

    log_info "Completed compile_per_protein_summary for ${PROTEIN}"
done

log_info "Done"
