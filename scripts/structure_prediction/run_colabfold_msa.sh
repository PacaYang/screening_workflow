#!/bin/bash
# Submit one ColabFold public-server MSA job per protein.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
PROTEINS="${MASTER_PROTEINS:-JAK1JH1}"
DRY_RUN=0

if [ "${MASTER_COLABFOLD_ENV+x}" = x ]; then
    # An explicitly empty value runs an executable that already has its own
    # environment, such as a full path supplied via MASTER_COLABFOLD_BATCH.
    COLABFOLD_ENV="$MASTER_COLABFOLD_ENV"
else
    COLABFOLD_ENV="colabfold"
fi
COLABFOLD_BATCH="${MASTER_COLABFOLD_BATCH:-colabfold_batch}"
PARTITION="${MASTER_COLABFOLD_PARTITION:-g24}"
CPUS="${MASTER_COLABFOLD_CPUS:-4}"
MEMORY="${MASTER_COLABFOLD_MEMORY:-8G}"
TIME_LIMIT="${MASTER_COLABFOLD_TIME_LIMIT:-08:00:00}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --task-root) TASK_ROOT="$2"; shift 2 ;;
        --proteins) PROTEINS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

SEQUENCES_CSV="${TASK_ROOT}/Input/sequences.csv"
[ -f "$SEQUENCES_CSV" ] || { echo "Sequences CSV not found: $SEQUENCES_CSV" >&2; exit 1; }

for protein in $PROTEINS; do
    msa_dir="${TASK_ROOT}/${protein}/fine_screening/MSA"
    log_dir="${msa_dir}/logs"
    job_ids="${msa_dir}/job_ids.txt"
    mkdir -p "$log_dir"
    : > "$job_ids"

    if python3 "${SCRIPT_DIR}/colabfold_msa.py" check \
        --sequences-csv "$SEQUENCES_CSV" --protein "$protein" \
        --msa-dir "$msa_dir" >/dev/null 2>&1; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: Reusing valid MSA cache for ${protein}"
        continue
    fi

    rm -f "${msa_dir}/msa.done"
    job_script="${log_dir}/slurm_colabfold_msa.sh"
    cat > "$job_script" <<EOF
#!/bin/bash
#SBATCH --job-name=msa_${protein}
#SBATCH --partition=${PARTITION}
#SBATCH --mem=${MEMORY}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${log_dir}/colabfold_msa_%j.out
#SBATCH --error=${log_dir}/colabfold_msa_%j.err

set -e
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
if [ -n "${COLABFOLD_ENV}" ]; then
    conda activate "${COLABFOLD_ENV}"
fi
if ! command -v "${COLABFOLD_BATCH}" >/dev/null 2>&1 && [ ! -x "${COLABFOLD_BATCH}" ]; then
    echo "ColabFold executable not found: ${COLABFOLD_BATCH}" >&2
    exit 1
fi

work_dir="${msa_dir}/work_\${SLURM_JOB_ID:-local}"
rm -rf "\$work_dir"
mkdir -p "\$work_dir"
python3 "${SCRIPT_DIR}/colabfold_msa.py" prepare \
    --sequences-csv "${SEQUENCES_CSV}" --protein "${protein}" \
    --work-dir "\$work_dir"

while IFS=$'\t' read -r representative fasta_path output_dir; do
    mkdir -p "\$output_dir"
    "${COLABFOLD_BATCH}" "\$fasta_path" "\$output_dir" --msa-only
done < "\$work_dir/search_plan.tsv"

colabfold_version=\$(python3 -c 'from importlib.metadata import version; print(version("colabfold"))' 2>/dev/null || echo unknown)
python3 "${SCRIPT_DIR}/colabfold_msa.py" finalize \
    --work-dir "\$work_dir" --msa-dir "${msa_dir}" \
    --colabfold-command "${COLABFOLD_BATCH} --msa-only" \
    --colabfold-version "\$colabfold_version"
rm -rf "\$work_dir"
EOF
    chmod +x "$job_script"

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "[DRY RUN] Would submit ColabFold MSA job for ${protein}: ${job_script}"
        continue
    fi
    job_id=$(sbatch --parsable "$job_script")
    [ -n "$job_id" ] || { echo "Failed to submit MSA job for ${protein}" >&2; exit 1; }
    echo "$job_id" >> "$job_ids"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: Submitted ColabFold MSA for ${protein}: ${job_id}"
done
