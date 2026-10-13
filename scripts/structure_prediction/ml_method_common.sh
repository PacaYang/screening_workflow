#!/bin/bash
# Shared execution helpers for NESSO, Protenix, and OpenFold3 wrappers.

set -e

ml_method_root() {
    local method="$1" upper
    upper=$(echo "$method" | tr '[:lower:]-' '[:upper:]_')
    local variable="MASTER_${upper}_ROOT"
    local value="${!variable:-}"
    printf '%s\n' "${value:-/home/yangl_pacagen_com/Applications/${method}}"
}

ml_method_env() {
    local method="$1" upper
    upper=$(echo "$method" | tr '[:lower:]-' '[:upper:]_')
    # An explicitly empty environment disables conda activation for wrappers
    # that are already running inside the desired environment.
    local variable="MASTER_${upper}_ENV"
    if [ "${!variable+x}" = x ]; then
        printf '%s\n' "${!variable}"
    else
        case "$method" in
            Nesso) printf '%s\n' 'nesso' ;;
            Protenix) printf '%s\n' 'protenix' ;;
            OpenFold3) printf '%s\n' 'openfold3' ;;
            *) printf '%s\n' "$method" ;;
        esac
    fi
}

ml_method_command() {
    local method="$1" upper
    upper=$(echo "$method" | tr '[:lower:]-' '[:upper:]_')
    local variable="MASTER_${upper}_CMD"
    local value="${!variable:-}"
    if [ -n "$value" ]; then
        printf '%s\n' "$value"
        return 0
    fi
    case "$method" in
        Nesso) printf '%s\n' 'nesso predict {input} --out_dir {output} --cache {weights} --accelerator gpu --override' ;;
        Protenix)
            if [ "${MASTER_PRECOMPUTE_MSA:-0}" -eq 1 ]; then
                printf '%s\n' 'protenix pred -i {input} -o {output} --use_msa true --sample 1 --step 20 --cycle 1'
            else
                printf '%s\n' 'protenix pred -i {input} -o {output} --use_msa false --sample 1 --step 20 --cycle 1'
            fi
            ;;
        OpenFold3) printf '%s\n' 'run_openfold predict --query_json {input} --output_dir {output} --inference_ckpt_path {weights}/home/.openfold3/of3-ob-2025-06-30-174k.pt --num_diffusion_samples 1 --num_model_seeds 1 --use_msa_server false' ;;
        *) return 1 ;;
    esac
}

ml_method_batches() {
    local method="$1" upper
    upper=$(echo "$method" | tr '[:lower:]-' '[:upper:]_')
    local variable="MASTER_${upper}_N_BATCHES"
    printf '%s\n' "${!variable:-40}"
}

ml_method_asset() {
    local method="$1" suffix="$2" default="$3" upper variable
    upper=$(echo "$method" | tr '[:lower:]-' '[:upper:]_')
    variable="MASTER_${upper}_${suffix}"
    printf '%s\n' "${!variable:-$default}"
}

ml_method_exe() {
    local method="$1" upper variable
    upper=$(echo "$method" | tr '[:lower:]-' '[:upper:]_')
    variable="MASTER_${upper}_EXE"
    printf '%s\n' "${!variable:-}"
}

ml_method_resource() {
    local method="$1" key="$2" default="$3" upper
    upper=$(echo "$method" | tr '[:lower:]-' '[:upper:]_')
    local variable="MASTER_${upper}_${key}"
    printf '%s\n' "${!variable:-$default}"
}

ml_method_run_command() {
    local method="$1" input="$2" output="$3" root command weights database exe
    root=$(ml_method_root "$method")
    exe=$(ml_method_exe "$method")
    weights=$(ml_method_asset "$method" WEIGHTS "/home/yangl_pacagen_com/Applications/model_weights/${method}")
    database=$(ml_method_asset "$method" DB "/home/yangl_pacagen_com/Applications/model_weights/${method}_db")
    command=$(ml_method_command "$method")
    command=${command//\{root\}/$root}
    command=${command//\{input\}/$input}
    command=${command//\{output\}/$output}
    command=${command//\{method\}/$method}
    command=${command//\{weights\}/$weights}
    command=${command//\{database\}/$database}
    command=${command//\{exe\}/$exe}
    echo "Running ${method}: ${command}"
    eval "$command"
}

ml_method_activate() {
    local method="$1" env weights
    env=$(ml_method_env "$method")
    if [ -n "$env" ]; then
        if [ -n "${CONDA_EXE:-}" ]; then
            # shellcheck disable=SC1091
            source "${CONDA_EXE%/bin/conda}/etc/profile.d/conda.sh" 2>/dev/null || true
        elif [ -f /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh ]; then
            # shellcheck disable=SC1091
            source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
        fi
        command -v conda >/dev/null 2>&1 || return 0
        conda activate "$env"
    fi
    weights=$(ml_method_asset "$method" WEIGHTS "/home/yangl_pacagen_com/Applications/model_weights/${method}")
    mkdir -p "$weights" "$weights/torch_extensions" "$weights/matplotlib"
    export TORCH_EXTENSIONS_DIR="$weights/torch_extensions"
    export MPLCONFIGDIR="$weights/matplotlib"
    case "$method" in
        Nesso) export NESSO_CACHE="$weights" ;;
        Protenix) export PROTENIX_ROOT_DIR="$weights/home" ;;
        OpenFold3) export OPENFOLD_CACHE="$weights/home/.openfold3" ;;
    esac
}
