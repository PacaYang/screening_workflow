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

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"
EXE="/home/yangl_pacagen_com/screening_workflow/scripts/input/gen_af3_json_template.py"
CONDA_ENV="HMSA_test"
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

log_info "Starting write_af3_template_input (n_ligands=${N_LIGANDS})"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    OUTDIR="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/input"
    TOKEN="${OUTDIR}/af3_input.done"
    SELECTED="${TASK_ROOT}/${PROTEIN}/initial_screening/selected.csv"
    # AF3 templates must be single-chain: <name>_<CHAIN>.cif lives in this dir.
    TEMPLATE_DIR="${TASK_ROOT}/Input/protein_file/${PROTEIN}"

    if [ -f "$TOKEN" ]; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    if [ ! -f "$SELECTED" ]; then
        log_error "selected.csv not found for ${PROTEIN}: $SELECTED"
        continue
    fi

    if [ ! -d "$TEMPLATE_DIR" ]; then
        log_error "Template dir not found for ${PROTEIN}: $TEMPLATE_DIR"
        continue
    fi

    conda activate "$CONDA_ENV"

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
