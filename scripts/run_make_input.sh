#!/bin/bash
#
# Standalone script for Snakemake checkpoint: make_input_csv
# Splits compound/sequence inputs into per-protein chunked CSVs.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"
SMILES_CSV="${TASK_ROOT}/Input/compounds_smiles.csv"
EXE="/home/yangl_pacagen_com/screening_workflow/scripts/make_input_csv.py"

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

log_info "Starting make_input_csv"

if [ ! -f "$SEQS_CSV" ]; then
    log_error "Sequences CSV not found: $SEQS_CSV"
    exit 1
fi

if [ ! -f "$SMILES_CSV" ]; then
    log_error "SMILES CSV not found: $SMILES_CSV"
    exit 1
fi

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    OUTDIR="${TASK_ROOT}/${PROTEIN}/initial_screening/inputs"
    TOKEN="${OUTDIR}/finish.token"

    if [ -f "$TOKEN" ]; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    conda activate HMSA_test

    mkdir -p "$OUTDIR"
    python "$EXE" \
        --sequences "$SEQS_CSV" \
        --smiles "$SMILES_CSV" \
        --protein "$PROTEIN" \
        --outdir "$OUTDIR"
    touch "$TOKEN"

    log_info "Completed make_input_csv for ${PROTEIN}"
done

log_info "Done"
