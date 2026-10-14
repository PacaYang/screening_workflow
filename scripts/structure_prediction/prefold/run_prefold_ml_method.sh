#!/bin/bash
# Generic protein-only prefold submitter for NESSO, Protenix, and OpenFold3.

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../ml_method_common.sh"
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

submit_prefold() {
    local protein="$1" method_dir="${TASK_ROOT}/${protein}/fine_screening/${METHOD}"
    local input_dir="${method_dir}/prefold" output_dir="${method_dir}/prefold/output"
    local log_dir="${method_dir}/prefold/logs" token="${method_dir}/prefold/prefold.done"
    local script="${log_dir}/slurm_prefold.sh" input_suffix="json" input
    [ "$METHOD" = "Nesso" ] && input_suffix="yaml"
    [ "$METHOD" = "Protenix" ] && input_suffix="protenix.json"
    [ "$METHOD" = "OpenFold3" ] && input_suffix="openfold.json"
    input="${input_dir}/${protein}.${input_suffix}"
    mkdir -p "$output_dir" "$log_dir"
    [ -f "$token" ] && return 0
    [ -f "$input" ] || { echo "Prefold input missing: $input" >&2; return 1; }
    local time_limit memory cpus partition
    time_limit=$(ml_method_resource "$METHOD" PREFOLD_TIME_LIMIT 48:00:00)
    memory=$(ml_method_resource "$METHOD" PREFOLD_MEMORY 64G)
    cpus=$(ml_method_resource "$METHOD" PREFOLD_CPUS 16)
    partition=$(ml_method_resource "$METHOD" PREFOLD_PARTITION g232)
    cat > "$script" <<EOF
#!/bin/bash
#SBATCH --job-name=${METHOD,,}_prefold_${protein}
#SBATCH --time=${time_limit}
#SBATCH --mem=${memory}
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --partition=${partition}
#SBATCH --output=${log_dir}/slurm_prefold_%j.out
#SBATCH --error=${log_dir}/slurm_prefold_%j.err
set -e
source "${SCRIPT_DIR}/../ml_method_common.sh"
ml_method_activate "${METHOD}"
mkdir -p "${output_dir}"
ml_method_run_command "${METHOD}" "${input}" "${output_dir}"
find "${output_dir}" -type f -print -quit | grep -q .
touch "${token}"
EOF
    chmod +x "$script"
    sbatch --parsable "$script" >> "${method_dir}/prefold/job_ids.txt"
}

for protein in $PROTEINS; do
    submit_prefold "$protein"
done
