#!/bin/bash
# Generic SLURM batch runner for NESSO, Protenix, and OpenFold3.

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/ml_method_common.sh"
METHOD=""
TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
PROTEINS="${MASTER_PROTEINS:-JAK1JH1}"
N_BATCHES=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --method) METHOD="$2"; shift 2 ;;
        --task-root) TASK_ROOT="$2"; shift 2 ;;
        --proteins) PROTEINS="$2"; shift 2 ;;
        --n-batches) N_BATCHES="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[ -n "$METHOD" ] || { echo "--method is required" >&2; exit 2; }
N_BATCHES="${N_BATCHES:-$(ml_method_batches "$METHOD")}";
TIME_LIMIT=$(ml_method_resource "$METHOD" TIME_LIMIT 48:00:00)
MEMORY=$(ml_method_resource "$METHOD" MEMORY 15G)
CPUS=$(ml_method_resource "$METHOD" CPUS 2)
PARTITION=$(ml_method_resource "$METHOD" PARTITION g24)

submit_batch() {
    local protein="$1" batch_id="$2"
    local method_dir="${TASK_ROOT}/${protein}/fine_screening/${METHOD}"
    local input_dir="${method_dir}/input" output_dir="${method_dir}/output" log_dir="${method_dir}/logs"
    local selected="${TASK_ROOT}/${protein}/initial_screening/selected.csv"
    local script="${log_dir}/slurm_batch_${batch_id}.sh"
    mkdir -p "$input_dir" "$output_dir/token" "$log_dir"
    cat > "$script" <<EOF
#!/bin/bash
#SBATCH --job-name=${METHOD,,}_${protein}_b${batch_id}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --mem=${MEMORY}
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --partition=${PARTITION}
#SBATCH --output=${log_dir}/slurm_batch_${batch_id}_%j.out
#SBATCH --error=${log_dir}/slurm_batch_${batch_id}_%j.err
set -e
source "${SCRIPT_DIR}/ml_method_common.sh"
ml_method_activate "${METHOD}"
source "${SCRIPT_DIR}/../pipeline/lib/job_result_ledger.sh"
selected="${selected}"
input_dir="${input_dir}"
output_dir="${output_dir}"
local_root="/tmp/${METHOD,,}_batch_${protein}_${batch_id}_\${SLURM_JOB_ID}"
mkdir -p "\$local_root"
n_smiles=\$(tail -n +2 "\$selected" | wc -l)
group_size=\${MASTER_N_LIGANDS:-1}
[ "\$group_size" -lt 1 ] && group_size=1
n_inputs=\$(( (n_smiles + group_size - 1) / group_size ))
batch_size=\$(( (n_inputs + ${N_BATCHES} - 1) / ${N_BATCHES} ))
start_idx=\$(( ${batch_id} * batch_size ))
end_idx=\$(( start_idx + batch_size ))
[ "\$end_idx" -gt "\$n_inputs" ] && end_idx="\$n_inputs"
success=0
failed=0
for idx in \$(seq "\$start_idx" \$((end_idx - 1))); do
    input_suffix="json"
    [ "${METHOD}" = "Nesso" ] && input_suffix="yaml"
    [ "${METHOD}" = "Protenix" ] && input_suffix="protenix.json"
    [ "${METHOD}" = "OpenFold3" ] && input_suffix="openfold.json"
    input="\${input_dir}/${protein}_\${idx}.\${input_suffix}"
    job_out="\${local_root}/job_\${idx}"
    mkdir -p "\$job_out"
    if [ ! -f "\$input" ]; then failed=\$((failed + 1)); continue; fi
    if ml_method_run_command "${METHOD}" "\$input" "\$job_out" && find "\$job_out" -type f -print -quit | grep -q .; then
        mkdir -p "\${output_dir}/job_\${idx}"
        cp -a "\$job_out/." "\${output_dir}/job_\${idx}/"
        primary_idx=\$((idx * group_size))
        python3 - "\$selected" "\$primary_idx" "\${output_dir}/job_\${idx}/smiles.txt" <<'PY'
import pandas as pd
import sys
df = pd.read_csv(sys.argv[1])
with open(sys.argv[3], "w", encoding="utf-8") as handle:
    handle.write(str(df.iloc[int(sys.argv[2])]["SMILES"]))
PY
        printf '%s\n' "${METHOD}" > "\${output_dir}/job_\${idx}/method.txt"
        touch "\${output_dir}/token/\${idx}.done"
        success=\$((success + 1))
    else
        failed=\$((failed + 1))
    fi
done
expected=\$((end_idx - start_idx))
status=failed
[ "\$expected" -gt 0 ] && [ "\$success" -eq "\$expected" ] && status=success
append_job_result "${output_dir}/job_results.tsv" "\${SLURM_JOB_ID:-unknown}" "${protein}" "batch_${batch_id}" "\$status" "\$success" "\$expected"
rm -rf "\$local_root"
[ "\$status" = success ]
EOF
    chmod +x "$script"
    sbatch --parsable "$script" >> "${output_dir}/job_ids.txt"
}

for protein in $PROTEINS; do
    selected="${TASK_ROOT}/${protein}/initial_screening/selected.csv"
    [ -f "$selected" ] || { echo "Selected CSV not found: $selected" >&2; continue; }
    output_dir="${TASK_ROOT}/${protein}/fine_screening/${METHOD}/output"
    mkdir -p "$output_dir"
    : > "${output_dir}/job_ids.txt"
    : > "${output_dir}/job_results.tsv"
    count=$(tail -n +2 "$selected" | wc -l)
    [ "$count" -gt 0 ] || continue
    group_size="${MASTER_N_LIGANDS:-1}"
    [ "$group_size" -lt 1 ] && group_size=1
    input_count=$(( (count + group_size - 1) / group_size ))
    batch_size=$(( (input_count + N_BATCHES - 1) / N_BATCHES ))
    for batch_id in $(seq 0 $((N_BATCHES - 1))); do
        [ $((batch_id * batch_size)) -lt "$input_count" ] || continue
        submit_batch "$protein" "$batch_id"
    done
done
