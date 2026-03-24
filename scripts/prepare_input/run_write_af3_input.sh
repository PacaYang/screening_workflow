#!/bin/bash
#
# Standalone script for Snakemake checkpoint: write_af3_input
# Generates AF3 JSON inputs with ligands for fine screening.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="${SCRIPT_DIR}/gen_af3_json_with_cmpds.py"

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

get_selected_csv() {
    local protein=$1
    local rel_path="${MASTER_SELECTED_REL_PATH:-initial_screening/selected.csv}"
    echo "${TASK_ROOT}/${protein}/${rel_path}"
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
    PREFOLD_JSON="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/prefold/${PROTEIN_LOWER}/${PROTEIN_LOWER}_data.json"
    SELECTED="$(get_selected_csv "${PROTEIN}")"
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
        log_error "Selected compounds file not found for ${PROTEIN}: $SELECTED"
        continue
    fi

    conda activate general

    echo "$PREFOLD_JSON"
    mkdir -p "$OUTDIR"
    python "$EXE" \
        --output-dir "$OUTDIR" \
        --input-json "$PREFOLD_JSON" \
        --smiles-file "$SELECTED"
    touch "$TOKEN"

    log_info "Completed write_af3_input for ${PROTEIN}"
done

log_info "Done"
