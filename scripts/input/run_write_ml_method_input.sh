#!/bin/bash
# Generic complex/template input writer for NESSO, Protenix, and OpenFold3.

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
PROTEINS="${MASTER_PROTEINS:-JAK1JH1}"
METHOD=""
MODE="complex"
N_LIGANDS="${MASTER_N_LIGANDS:-1}"
PRECOMPUTE_MSA="${MASTER_PRECOMPUTE_MSA:-0}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --method) METHOD="$2"; shift 2 ;;
        --task-root) TASK_ROOT="$2"; shift 2 ;;
        --proteins) PROTEINS="$2"; shift 2 ;;
        --mode) MODE="$2"; shift 2 ;;
        --n-ligands) N_LIGANDS="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[ -n "$METHOD" ] || { echo "--method is required" >&2; exit 2; }

for protein in $PROTEINS; do
    method_dir="${TASK_ROOT}/${protein}/fine_screening/${METHOD}"
    out_dir="${method_dir}/input"
    token="${method_dir}/${METHOD,,}_input.done"
    mode_file="${method_dir}/${METHOD,,}_input.mode"
    selected="${TASK_ROOT}/${protein}/initial_screening/selected.csv"
    [ -f "$selected" ] || { echo "Selected CSV not found: $selected" >&2; continue; }
    expected_mode="$MODE"
    msa_args=()
    if [ "$PRECOMPUTE_MSA" -eq 1 ] && { [ "$METHOD" = Protenix ] || [ "$METHOD" = OpenFold3 ]; }; then
        msa_manifest="${TASK_ROOT}/${protein}/fine_screening/MSA/manifest.json"
        msa_token="${TASK_ROOT}/${protein}/fine_screening/MSA/msa.done"
        [ -f "$msa_manifest" ] && [ -f "$msa_token" ] || {
            echo "Shared MSA is incomplete for ${protein}: ${msa_manifest}" >&2
            exit 1
        }
        expected_mode="precomputed_msa:$(sha256sum "$msa_manifest" | awk '{print $1}')"
        msa_args+=(--msa-manifest "$msa_manifest")
    fi
    current_mode="complex"
    [ -f "$mode_file" ] && current_mode=$(cat "$mode_file")
    if [ -f "$token" ] && [ "$current_mode" = "$expected_mode" ]; then
        continue
    fi
    mkdir -p "$out_dir"
    template_args=()
    if [ "$MODE" = template ]; then
        template_args+=(--template-dir "${TASK_ROOT}/Input/protein_file/${protein}")
    fi
    python3 "${SCRIPT_DIR}/gen_ml_method_input.py" \
        --method "$METHOD" --mode "$MODE" --output-dir "$out_dir" \
        --input-csv "${TASK_ROOT}/Input/sequences.csv" --protein-name "$protein" \
        --smiles-file "$selected" --n-ligands "$N_LIGANDS" \
        "${template_args[@]}" "${msa_args[@]}"
    printf '%s\n' "$expected_mode" > "$mode_file"
    touch "$token"
done
