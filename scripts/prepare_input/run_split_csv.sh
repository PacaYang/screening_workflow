#!/bin/bash
#
# Standalone script for Snakemake checkpoint: split_csv
# Splits selected compounds into chunked CSVs for Vina docking.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="${SCRIPT_DIR}/split_csv.py"

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

log_info "Starting split_csv"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    SELECTED="$(get_selected_csv "${PROTEIN}")"
    PDB="${TASK_ROOT}/Input/protein_file/${PROTEIN}/${PROTEIN}.pdb"
    OUTDIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/input"

    # Skip if output dir already has CSV files
    if ls "${OUTDIR}"/*.csv &>/dev/null; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    if [ ! -f "$SELECTED" ]; then
        log_error "Selected compounds file not found for ${PROTEIN}: $SELECTED"
        continue
    fi

    if [ ! -f "$PDB" ]; then
        log_error "PDB file not found for ${PROTEIN}: $PDB"
        continue
    fi

    conda activate general

    mkdir -p "$OUTDIR"
    python "$EXE" \
        --protein-name="$PROTEIN" \
        --pdb-file="$PDB" \
        --smiles-file="$SELECTED" \
        --chunk-size=100 \
        --output-dir="$OUTDIR"

    log_info "Completed split_csv for ${PROTEIN}"
done

log_info "Done"
