#!/bin/bash
#
# Standalone script for Snakemake checkpoint: write_boltz2_input
# Generates Boltz2 YAML inputs with ligands for fine screening.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"
EXE="${SCRIPT_DIR}/gen_boltz_yaml.py"
PRECOMPUTE_MSA="${MASTER_PRECOMPUTE_MSA:-0}"
N_LIGANDS="${MASTER_N_LIGANDS:-1}"

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

log_info "Starting write_boltz2_input"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    OUTDIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/input"
    TOKEN="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/boltz_input.done"
    MODE_FILE="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/boltz_input.mode"
    MSA="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/prefold/boltz_results_${PROTEIN}/msa/${PROTEIN}_0.csv"
    SELECTED="${TASK_ROOT}/${PROTEIN}/initial_screening/selected.csv"
    CONFIDENCE="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/prefold/boltz_results_${PROTEIN}/predictions/${PROTEIN}/confidence_${PROTEIN}_model_0.json"

    EXPECTED_MODE="legacy"
    if [ "$PRECOMPUTE_MSA" -eq 1 ]; then
        MANIFEST="${TASK_ROOT}/${PROTEIN}/fine_screening/MSA/manifest.json"
        MSA_TOKEN="${TASK_ROOT}/${PROTEIN}/fine_screening/MSA/msa.done"
        [ -f "$MANIFEST" ] && [ -f "$MSA_TOKEN" ] || {
            log_error "Shared MSA is incomplete for ${PROTEIN}: ${MANIFEST}"
            exit 1
        }
        EXPECTED_MODE="precomputed_msa:$(sha256sum "$MANIFEST" | awk '{print $1}')"
    fi
    CURRENT_MODE="legacy"
    [ -f "$MODE_FILE" ] && CURRENT_MODE=$(cat "$MODE_FILE")
    if [ -f "$TOKEN" ] && [ "$CURRENT_MODE" = "$EXPECTED_MODE" ]; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    if [ "$PRECOMPUTE_MSA" -eq 1 ]; then
        [ -f "$SELECTED" ] || { log_error "selected.csv not found: $SELECTED"; exit 1; }
        conda activate boltz_test
        mkdir -p "$OUTDIR"
        python "$EXE" \
            --output "$OUTDIR" \
            --msa-manifest "$MANIFEST" \
            --smiles-path "$SELECTED" \
            --protein-name "$PROTEIN" \
            --protein-file "$SEQS_CSV" \
            --name-col protein_name \
            --n-ligands "$N_LIGANDS"
        printf '%s\n' "$EXPECTED_MODE" > "$MODE_FILE"
        touch "$TOKEN"
        log_info "Completed precomputed-MSA Boltz2 inputs for ${PROTEIN}"
        continue
    fi

    if [ ! -f "$CONFIDENCE" ]; then
        log_error "Boltz2 prefold confidence JSON not found for ${PROTEIN}: $CONFIDENCE"
        continue
    fi

    if [ ! -f "$MSA" ]; then
        log_error "Boltz2 prefold MSA not found for ${PROTEIN}: $MSA"
        continue
    fi

    if [ ! -f "$SELECTED" ]; then
        log_error "selected.csv not found for ${PROTEIN}: $SELECTED"
        continue
    fi

    conda activate boltz_test

    mkdir -p "$OUTDIR"
    python "$EXE" \
        --output "$OUTDIR" \
        --msa "$MSA" \
        --smiles-path "$SELECTED" \
        --protein-name "$PROTEIN" \
        --protein-file "$SEQS_CSV" \
        --name-col protein_name
    printf '%s\n' "$EXPECTED_MODE" > "$MODE_FILE"
    touch "$TOKEN"

    log_info "Completed write_boltz2_input for ${PROTEIN}"
done

log_info "Done"
