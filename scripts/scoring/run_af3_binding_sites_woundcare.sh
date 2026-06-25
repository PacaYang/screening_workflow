#!/bin/bash
#SBATCH --job-name=af3_binding_sites_woundcare
#SBATCH --output=/shared/cuteness_woundcare/logs/af3_binding_sites_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --priority=TOP

set -e

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="/shared/cuteness_woundcare"
SCRIPTS_ROOT="/home/ubuntu/screening_workflow/scripts"
GENERAL_PYTHON="/home/ubuntu/miniconda3/envs/general/bin/python"

EXCLUDE_REGEX='^(JAK1|JAK2|JAK3|TYK2|JAK_input|Input|logs)$'

mkdir -p "${TASK_ROOT}/logs"

PROTEINS=()
for entry in "${TASK_ROOT}"/*; do
    [ -d "$entry" ] || continue
    name=$(basename "$entry")
    if [[ "$name" =~ $EXCLUDE_REGEX ]]; then
        continue
    fi
    if [ -d "${entry}/fine_screening/AF3/output" ]; then
        PROTEINS+=("$name")
    fi
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found ${#PROTEINS[@]} proteins to process"

for protein in "${PROTEINS[@]}"; do
    af3_dir="${TASK_ROOT}/${protein}/fine_screening/AF3"
    af3_out="${af3_dir}/output"
    summary_csv="${af3_dir}/summary.csv"

    if [ -f "$summary_csv" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skip ${protein}: summary.csv exists"
        continue
    fi

    archives=("${af3_out}"/batch_*.tar.gz)
    if [ ! -e "${archives[0]}" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skip ${protein}: no batch archives"
        continue
    fi

    scratch=$(mktemp -d -t "af3_${protein}_XXXXXX")
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${protein}: extracting to ${scratch}"
    for tgz in "${archives[@]}"; do
        tar -xzf "$tgz" -C "$scratch"
    done

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${protein}: scoring"
    "$GENERAL_PYTHON" "${SCRIPTS_ROOT}/scoring/af3_scores.py" \
        --af3-results-folder "$scratch" \
        --output-dir "$af3_dir" \
        --binding-sites \
        --distance-threshold 10.0

    rm -rf "$scratch"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done: ${summary_csv}"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All proteins complete."
