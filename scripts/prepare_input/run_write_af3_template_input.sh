#!/bin/bash
#
# Template-mode input writer for AF3 (single stage, no prefold).
# Generates AF3 JSON inputs with per-chain templates + empty MSA, multi-chain
# and multi-ligand, reading sequences directly from Input/sequences.csv.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"
EXE="${SCRIPT_DIR}/gen_af3_json_template.py"
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

get_selected_csv() {
    local protein=$1
    local rel_path="${MASTER_SELECTED_REL_PATH:-initial_screening/selected.csv}"
    echo "${TASK_ROOT}/${protein}/${rel_path}"
}

# ============================================================================
# Main
# ============================================================================

log_info "Starting write_af3_template_input (n_ligands=${N_LIGANDS})"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    OUTDIR="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/input"
    TOKEN="${OUTDIR}/af3_input.done"
    SELECTED="$(get_selected_csv "${PROTEIN}")"
    # AF3 templates must be single-chain: <name>_<CHAIN>.cif lives in this dir.
    TEMPLATE_DIR="${TASK_ROOT}/Input/protein_file/${PROTEIN}"

    if [ -f "$TOKEN" ]; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    if [ ! -f "$SELECTED" ]; then
        log_error "Selected compounds file not found for ${PROTEIN}: $SELECTED"
        continue
    fi

    if [ ! -d "$TEMPLATE_DIR" ]; then
        log_error "Template dir not found for ${PROTEIN}: $TEMPLATE_DIR"
        continue
    fi

    conda activate general

    mkdir -p "$OUTDIR"
    python "$EXE" \
        --output-dir "$OUTDIR" \
        --input-csv "$SEQS_CSV" \
        --protein-name "$PROTEIN" \
        --smiles-file "$SELECTED" \
        --template-dir "$TEMPLATE_DIR" \
        --n-ligands "$N_LIGANDS"
    touch "$TOKEN"

    log_info "Completed write_af3_template_input for ${PROTEIN}"
done

log_info "Done"
