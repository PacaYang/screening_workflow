#!/bin/bash
#
# Standalone script for Snakemake checkpoint: write_boltz2_input
# Generates Boltz2 YAML inputs with ligands for fine screening.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"
EXE="${SCRIPT_DIR}/gen_boltz_yaml.py"

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

log_info "Starting write_boltz2_input"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    OUTDIR="/home/ubuntu/${PROTEIN}/boltz2_tmp/input"
    TOKEN="/home/ubuntu/${PROTEIN}/boltz2_tmp/boltz_input.done"
    MSA="/home/ubuntu/${PROTEIN}/boltz2_tmp/boltz_results_${PROTEIN}/msa/${PROTEIN}_0.csv"
    SELECTED="$(get_selected_csv "${PROTEIN}")"
    CONFIDENCE="/home/ubuntu/${PROTEIN}/boltz2_tmp/boltz_results_${PROTEIN}/predictions/${PROTEIN}/confidence_${PROTEIN}_model_0.json"

    if [ -f "$TOKEN" ]; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    if [ ! -f "$CONFIDENCE" ]; then
        log_error "Boltz2 prefold confidence JSON not found for ${PROTEIN}: $CONFIDENCE"
        continue
    fi

    if [ ! -f "$SELECTED" ]; then
        log_error "Selected compounds file not found for ${PROTEIN}: $SELECTED"
        continue
    fi

    conda activate boltz

    mkdir -p "$OUTDIR"
    python "$EXE" \
        --output "$OUTDIR" \
        --msa "$MSA" \
        --smiles-path "$SELECTED" \
        --protein-name "$PROTEIN" \
        --protein-file "$SEQS_CSV"
    touch "$TOKEN"

    log_info "Completed write_boltz2_input for ${PROTEIN}"
done

log_info "Done"
