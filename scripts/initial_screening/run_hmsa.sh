#!/bin/bash
#
# Standalone script for Snakemake rule: run_hmsa
# Submits SLURM jobs for HMSA initial screening (one per input chunk).
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
EXE="/home/yangl_pacagen_com/screening_workflow/scripts/initial_screening/HMSA_predict.py"
MODEL="${MASTER_HMSA_MODEL:-/home/yangl_pacagen_com/Applications/model_weights/HMSA/model.pt}"
CONDA_ENV="HMSA_test"

# SLURM configuration
TIME_LIMIT="48:00:00"
MEMORY="15G"
CPUS_PER_TASK=2
PARTITION="g24"

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

submit_job() {
    local protein=$1
    local idx=$2

    local INIT_DIR="${TASK_ROOT}/${protein}/initial_screening"
    local INPUT_CSV="${INIT_DIR}/inputs/input_${idx}.csv"
    local OUTPUT_DIR="${INIT_DIR}/HMSA"
    local OUTPUT_CSV="${OUTPUT_DIR}/prediction_${idx}.csv"

    if [ -f "$OUTPUT_CSV" ]; then
        log_info "prediction_${idx}.csv already exists for ${protein}, skipping"
        return 0
    fi

    mkdir -p "$OUTPUT_DIR"

    local JOB_SCRIPT="${OUTPUT_DIR}/slurm_${idx}.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=hmsa_${protein}_${idx}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --mem=${MEMORY}
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --partition=${PARTITION}
#SBATCH --output=${OUTPUT_DIR}/slurm_${idx}_%j.out
#SBATCH --error=${OUTPUT_DIR}/slurm_${idx}_%j.err

set -e

echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

mkdir -p "${OUTPUT_DIR}"

export PYTHONPATH="/home/yangl_pacagen_com/Applications/HMSA-DTI:\${PYTHONPATH}"
python "${EXE}" \\
    --test_path "${INPUT_CSV}" \\
    --preds_path "${OUTPUT_CSV}" \\
    --checkpoint_paths "${MODEL}" \\
    --smiles_columns "SMILES"

echo "Job completed at: \$(date)"
EOF

    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted input_${idx} for ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${OUTPUT_DIR}/job_ids.txt"
        return 0
    else
        log_error "Failed to submit input_${idx} for ${protein}"
        return 1
    fi
}

# ============================================================================
# Main
# ============================================================================

log_info "Starting HMSA batch job submission"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    INPUT_DIR="${TASK_ROOT}/${PROTEIN}/initial_screening/inputs"
    TOKEN="${INPUT_DIR}/finish.token"

    if [ ! -f "$TOKEN" ]; then
        log_error "finish.token not found for ${PROTEIN}: ${TOKEN}"
        log_error "Please run run_make_input.sh first"
        continue
    fi

    INPUT_FILES=$(ls "${INPUT_DIR}"/input_*.csv 2>/dev/null | sort)

    if [ -z "$INPUT_FILES" ]; then
        log_error "No input CSV files found in ${INPUT_DIR}"
        continue
    fi

    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/initial_screening/HMSA"
    mkdir -p "$OUTPUT_DIR"
    > "${OUTPUT_DIR}/job_ids.txt"

    SUBMITTED=0
    for INPUT_FILE in $INPUT_FILES; do
        IDX=$(basename "$INPUT_FILE" .csv | sed 's/input_//')
        if submit_job "$PROTEIN" "$IDX"; then
            SUBMITTED=$((SUBMITTED + 1))
        fi
    done

    log_info "Submitted ${SUBMITTED} HMSA jobs for ${PROTEIN}"
done

log_info "All HMSA jobs submitted"
log_info "Monitor jobs with: squeue -u \$USER"
