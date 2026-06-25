#!/bin/bash
#
# Template-mode input writer for Boltz2 (single stage, no prefold).
# Generates Boltz2 YAML inputs with template + empty MSA, multi-chain and
# multi-ligand, reading sequences directly from Input/sequences.csv.
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

log_info "Starting write_boltz2_template_input (n_ligands=${N_LIGANDS})"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    OUTDIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/input"
    TOKEN="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/boltz_input.done"
    SELECTED="$(get_selected_csv "${PROTEIN}")"
    # Multi-chain mmCIF template (Boltz accepts multiple chains in one file).
    TEMPLATE_CIF="${TASK_ROOT}/Input/protein_file/${PROTEIN}/${PROTEIN}.cif"

    if [ -f "$TOKEN" ]; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    if [ ! -f "$SELECTED" ]; then
        log_error "Selected compounds file not found for ${PROTEIN}: $SELECTED"
        continue
    fi

    if [ ! -f "$TEMPLATE_CIF" ]; then
        log_error "Template CIF not found for ${PROTEIN}: $TEMPLATE_CIF"
        continue
    fi

    conda activate boltz

    mkdir -p "$OUTDIR"
    python "$EXE" \
        --output "$OUTDIR" \
        --template-mode \
        --template-cif "$TEMPLATE_CIF" \
        --n-ligands "$N_LIGANDS" \
        --smiles-path "$SELECTED" \
        --protein-name "$PROTEIN" \
        --protein-file "$SEQS_CSV" \
        --name-col name
    touch "$TOKEN"

    log_info "Completed write_boltz2_template_input for ${PROTEIN}"
done

log_info "Done"
