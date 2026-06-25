#!/bin/bash
#
# Prepare AutoDock Vina receptor PDBQT using obabel.
#

set -euo pipefail

TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
PROTEIN=""
OUTPUT_FILE=""
CONDA_ENV=""
FORCE=0

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

usage() {
    cat <<USAGE
Usage: $0 --protein NAME [--task-root DIR] [--output FILE] [--conda-env ENV] [--force]

Options:
  --protein NAME     Protein name (required)
  --task-root DIR    Task root (default: MASTER_TASK_ROOT or /home/ubuntu/snake_test)
  --output FILE      Output receptor PDBQT path
  --conda-env ENV    Activate this conda env before running obabel
  --force            Rebuild receptor even if output exists
  --help             Show help
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --protein)
            PROTEIN="$2"
            shift 2
            ;;
        --task-root)
            TASK_ROOT="$2"
            shift 2
            ;;
        --output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [ -z "$PROTEIN" ]; then
    log_error "Missing required argument: --protein"
    usage
    exit 1
fi

INPUT_PDB="${TASK_ROOT}/Input/protein_file/${PROTEIN}/${PROTEIN}.pdb"

if [ -z "$OUTPUT_FILE" ]; then
    OUTPUT_FILE="${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/receptor/${PROTEIN}.pdbqt"
fi

if [ ! -f "$INPUT_PDB" ]; then
    log_error "Protein PDB not found: ${INPUT_PDB}"
    exit 1
fi

if [ -n "$CONDA_ENV" ]; then
    source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
    conda activate "$CONDA_ENV"
fi

if ! command -v obabel >/dev/null 2>&1; then
    log_error "obabel is not available in PATH. Activate the correct environment and retry."
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_FILE")"

if [ -f "$OUTPUT_FILE" ] && [ "$FORCE" -eq 0 ]; then
    log_info "Reusing existing receptor PDBQT: ${OUTPUT_FILE}"
    echo "$OUTPUT_FILE"
    exit 0
fi

log_info "Preparing receptor with obabel: ${INPUT_PDB} -> ${OUTPUT_FILE}"
if ! obabel -ipdb "$INPUT_PDB" -opdbqt -O "$OUTPUT_FILE" -xr --addpolarh --partialcharge gasteiger; then
    log_error "obabel receptor preparation failed for ${PROTEIN}"
    exit 1
fi

log_info "Receptor preparation completed: ${OUTPUT_FILE}"
echo "$OUTPUT_FILE"
