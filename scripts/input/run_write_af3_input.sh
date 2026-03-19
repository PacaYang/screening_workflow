#!/bin/bash
#
# Standalone script for Snakemake checkpoint: write_af3_input
# Generates AF3 JSON inputs with ligands for fine screening.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
EXE="/home/yangl_pacagen_com/screening_workflow/scripts/input/gen_af3_json_with_cmpds.py"

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

log_info "Starting write_af3_input"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    PROTEIN_LOWER=$(echo "$PROTEIN" | tr '[:upper:]' '[:lower:]')
    OUTDIR="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/input"
    TOKEN="${OUTDIR}/af3_input.done"
    PREFOLD_JSON="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/prefold/${PROTEIN}/${PROTEIN}_data.json"
    SELECTED="${TASK_ROOT}/${PROTEIN}/initial_screening/selected.csv"
    PREFOLD_TOKEN="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/prefold/prefold.done"

    if [ -f "$TOKEN" ]; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    if [ ! -f "$PREFOLD_TOKEN" ]; then
        log_error "prefold.done not found for ${PROTEIN}: $PREFOLD_TOKEN"
        continue
    fi

    if [ ! -f "$SELECTED" ]; then
        log_error "selected.csv not found for ${PROTEIN}: $SELECTED"
        continue
    fi

    conda activate HMSA_test

    echo "$PREFOLD_JSON"
    mkdir -p "$OUTDIR"
    python "$EXE" \
        --output-dir "$OUTDIR" \
        --input-json "$PREFOLD_JSON" \
        --smiles-file "$SELECTED" \
        --smiles-col "SMILES"
    touch "$TOKEN"

    log_info "Completed write_af3_input for ${PROTEIN}"
done

log_info "Done"
