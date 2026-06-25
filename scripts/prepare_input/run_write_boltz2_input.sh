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

    OUTDIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/input"
    TOKEN="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/boltz_input.done"
    BOLTZ_TMP="/home/ubuntu/${PROTEIN}/boltz2_tmp"
    MSA="${BOLTZ_TMP}/boltz_results_${PROTEIN}/msa/${PROTEIN}_0.csv"
    SELECTED="$(get_selected_csv "${PROTEIN}")"
    CONFIDENCE="${BOLTZ_TMP}/boltz_results_${PROTEIN}/predictions/${PROTEIN}/confidence_${PROTEIN}_model_0.json"

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

    # Use prefold YAML as template to preserve oligomeric state (e.g., trimers)
    PREFOLD_YAML="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/prefold/${PROTEIN}.yaml"

    if [ -f "$PREFOLD_YAML" ]; then
        log_info "Using prefold YAML template: $PREFOLD_YAML"
        python "$EXE" \
            --output "$OUTDIR" \
            --prefold-yaml "$PREFOLD_YAML" \
            --msa "$MSA" \
            --smiles-path "$SELECTED"
    else
        log_info "Prefold YAML not found, using legacy single-chain mode"
        python "$EXE" \
            --output "$OUTDIR" \
            --msa "$MSA" \
            --smiles-path "$SELECTED" \
            --protein-name "$PROTEIN" \
            --protein-file "$SEQS_CSV" \
            --name-col name
    fi
    touch "$TOKEN"

    log_info "Completed write_boltz2_input for ${PROTEIN}"
done

log_info "Done"
