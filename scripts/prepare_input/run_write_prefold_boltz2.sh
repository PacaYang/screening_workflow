#!/bin/bash
#
# Standalone script for Snakemake rule: write_prefold_boltz2
# Generates Boltz2 YAML input for protein-only prefold.
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

# ============================================================================
# Main
# ============================================================================

log_info "Starting write_prefold_boltz2"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    OUTDIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/prefold"
    OUTPUT_YAML="${OUTDIR}/${PROTEIN}.yaml"

    if [ -f "$OUTPUT_YAML" ]; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    conda activate boltz

    mkdir -p "$OUTDIR"
    python "$EXE" \
        --output "$OUTDIR" \
        --protein-name "$PROTEIN" \
        --protein-file "$SEQS_CSV" \
        --protein-only

    log_info "Completed write_prefold_boltz2 for ${PROTEIN}"
done

log_info "Done"
