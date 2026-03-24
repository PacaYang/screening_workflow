#!/bin/bash
#
# Standalone script for Snakemake rule: run_colddta
# Submits SLURM jobs for ColdDTA initial screening (one per input chunk).
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
EXE="/home/ubuntu/screening_workflow/algos/coldDTA/predict.py"
CHECKPOINT="/home/ubuntu/screening_workflow/algos/coldDTA/model/epoch1297test_loss0.1798.pt"
CONDA_ENV="cold"

# SLURM configuration
TIME_LIMIT="48:00:00"
MEMORY="15G"
CPUS_PER_TASK=4
CONSTRAINT="g5.xlarge"

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
    local OUTPUT_DIR="${INIT_DIR}/ColdDTA"
    local OUTPUT_CSV="${OUTPUT_DIR}/prediction_${idx}.csv"
    local FAILED_CSV="${OUTPUT_DIR}/failed_smiles_${idx}.csv"

    if [ -f "$OUTPUT_CSV" ]; then
        log_info "prediction_${idx}.csv already exists for ${protein}, skipping"
        return 0
    fi

    mkdir -p "$OUTPUT_DIR"

    local JOB_SCRIPT="${OUTPUT_DIR}/slurm_${idx}.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=colddta_${protein}_${idx}
#SBATCH --constraint=${CONSTRAINT}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --mem=${MEMORY}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --output=${OUTPUT_DIR}/slurm_${idx}_%j.out
#SBATCH --error=${OUTPUT_DIR}/slurm_${idx}_%j.err

set -e

echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

mkdir -p "${OUTPUT_DIR}"

python "${EXE}" \\
    --input "${INPUT_CSV}" \\
    --output "${OUTPUT_CSV}" \\
    --failed-smiles "${FAILED_CSV}" \\
    --batch-size 64 \\
    --target-length 1000 \\
    --checkpoint "${CHECKPOINT}"

echo "Job completed at: \$(date)"
EOF

    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted input_${idx} for ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${OUTPUT_DIR}/submitted_jobs.txt"
        return 0
    else
        log_error "Failed to submit input_${idx} for ${protein}"
        return 1
    fi
}

# ============================================================================
# Main
# ============================================================================

log_info "Starting ColdDTA batch job submission"

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

    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/initial_screening/ColdDTA"
    mkdir -p "$OUTPUT_DIR"
    > "${OUTPUT_DIR}/submitted_jobs.txt"

    SUBMITTED=0
    for INPUT_FILE in $INPUT_FILES; do
        IDX=$(basename "$INPUT_FILE" .csv | sed 's/input_//')
        if submit_job "$PROTEIN" "$IDX"; then
            SUBMITTED=$((SUBMITTED + 1))
        fi
    done

    log_info "Submitted ${SUBMITTED} ColdDTA jobs for ${PROTEIN}"
done

log_info "All ColdDTA jobs submitted"
log_info "Monitor jobs with: squeue -u \$USER"
