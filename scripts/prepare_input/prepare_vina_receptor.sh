#!/bin/bash
#
# Prepare AutoDock Vina receptor PDBQT using obabel.
# Requires ROOT/ENDROOT markers in the generated receptor file.
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

validate_root_block() {
    local pdbqt_file=$1
    local root_line
    local end_line
    local between_lines

    root_line=$(awk '/^[[:space:]]*ROOT[[:space:]]*$/ {print NR; exit}' "$pdbqt_file")
    end_line=$(awk '/^[[:space:]]*ENDROOT[[:space:]]*$/ {line=NR} END {if (line) print line}' "$pdbqt_file")

    if [ -z "${root_line:-}" ] || [ -z "${end_line:-}" ]; then
        log_error "Missing ROOT/ENDROOT markers in receptor PDBQT: ${pdbqt_file}"
        return 1
    fi

    if [ "$end_line" -le "$root_line" ]; then
        log_error "Invalid ROOT/ENDROOT ordering in receptor PDBQT: ${pdbqt_file}"
        return 1
    fi

    between_lines=$((end_line - root_line - 1))
    if [ "$between_lines" -le 0 ]; then
        log_error "No content found between ROOT and ENDROOT in receptor PDBQT: ${pdbqt_file}"
        return 1
    fi

    return 0
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
    if validate_root_block "$OUTPUT_FILE"; then
        log_info "Reusing existing receptor PDBQT: ${OUTPUT_FILE}"
        echo "$OUTPUT_FILE"
        exit 0
    fi
    log_info "Existing receptor PDBQT is invalid. Rebuilding: ${OUTPUT_FILE}"
fi

log_info "Preparing receptor with obabel: ${INPUT_PDB} -> ${OUTPUT_FILE}"
if ! obabel -ipdb "$INPUT_PDB" -opdbqt -O "$OUTPUT_FILE" -xr --addpolarh --partialcharge gasteiger; then
    log_error "obabel receptor preparation failed for ${PROTEIN}"
    exit 1
fi

if ! validate_root_block "$OUTPUT_FILE"; then
    exit 1
fi

log_info "Receptor preparation completed: ${OUTPUT_FILE}"
echo "$OUTPUT_FILE"
