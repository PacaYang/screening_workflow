#!/bin/bash
# Generic protein-only input writer for NESSO, Protenix, and OpenFold3.

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
PROTEINS="${MASTER_PROTEINS:-JAK1JH1}"
METHOD=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --method) METHOD="$2"; shift 2 ;;
        --task-root) TASK_ROOT="$2"; shift 2 ;;
        --proteins) PROTEINS="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[ -n "$METHOD" ] || { echo "--method is required" >&2; exit 2; }

for protein in $PROTEINS; do
    out_dir="${TASK_ROOT}/${protein}/fine_screening/${METHOD}/prefold"
    token="${out_dir}/prefold_input.done"
    [ -f "$token" ] && continue
    mkdir -p "$out_dir"
    python3 "${SCRIPT_DIR}/gen_ml_method_input.py" \
        --method "$METHOD" --mode prefold --output-dir "$out_dir" \
        --input-csv "${TASK_ROOT}/Input/sequences.csv" --protein-name "$protein"
    touch "$token"
done
