#!/bin/bash
#
# Standalone script for Snakemake rule: write_prefold_af3
# Generates AF3 JSON input for protein-only prefold.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"
EXE="/home/ubuntu/screening_workflow/scripts/gen_af3_json_protein.py"

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

log_info "Starting write_prefold_af3"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    OUTDIR="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/prefold"
    OUTPUT_JSON="${OUTDIR}/${PROTEIN}.json"

    if [ -f "$OUTPUT_JSON" ]; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    mkdir -p "$OUTDIR"
    python "$EXE" \
        --output-dir "$OUTDIR" \
        --protein-name "$PROTEIN" \
        --input-csv "$SEQS_CSV"

    log_info "Completed write_prefold_af3 for ${PROTEIN}"
done

log_info "Done"
