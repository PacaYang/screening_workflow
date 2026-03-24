#!/bin/bash
#
# Standalone script for Snakemake rule: prefold_af3
# Submits a SLURM job per protein to run AlphaFold3 protein-only folding.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
AF3_EXE="/home/ubuntu/Applications/alphafold3/run_alphafold.py"
AF3_WEIGHT_DIR="/shared/programs/af3_weights"
AF3_DB_DIR="/shared/programs/af3_data"
HMMER_DIR="/home/ubuntu/Applications/hmmer/bin"
CONDA_ENV="af3"

# SLURM configuration
TIME_LIMIT="48:00:00"
MEMORY="15G"
CPUS_PER_TASK=4
CONSTRAINT="g5.xlarge"

# ============================================================================
# Functions
# ============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*" >&2
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

submit_prefold_job() {
    local protein=$1
    local protein_lower=$(echo "$protein" | tr '[:upper:]' '[:lower:]')

    local PREFOLD_DIR="${TASK_ROOT}/${protein}/fine_screening/AF3/prefold"
    local INPUT_JSON="${PREFOLD_DIR}/${protein}.json"
    local TOKEN="${PREFOLD_DIR}/prefold.done"
    local LOG_DIR="${PREFOLD_DIR}/logs"

    if [ -f "$TOKEN" ]; then
        log_info "prefold.done already exists for ${protein}, skipping"
        return 0
    fi

    if [ ! -f "$INPUT_JSON" ]; then
        log_error "Input JSON not found for ${protein}: ${INPUT_JSON}"
        log_error "Please run run_write_prefold_af3.sh first"
        return 1
    fi

    mkdir -p "$LOG_DIR"

    local JOB_SCRIPT="${LOG_DIR}/slurm_prefold.sh"

    cat > "$JOB_SCRIPT" <<'EOFSCRIPT'
#!/bin/bash
#SBATCH --job-name=af3_prefold_PROTEIN_PLACEHOLDER
#SBATCH --constraint=CONSTRAINT_PLACEHOLDER
#SBATCH --time=TIME_LIMIT_PLACEHOLDER
#SBATCH --mem=MEMORY_PLACEHOLDER
#SBATCH --cpus-per-task=CPUS_PLACEHOLDER
#SBATCH --output=LOG_DIR_PLACEHOLDER/slurm_prefold_%j.out
#SBATCH --error=LOG_DIR_PLACEHOLDER/slurm_prefold_%j.err

set +eu

echo "Job started at: $(date)"
echo "Running on host: $(hostname)"
echo "Folding protein: PROTEIN_PLACEHOLDER"

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate CONDA_ENV_PLACEHOLDER

python -u AF3_EXE_PLACEHOLDER \
    --json_path INPUT_JSON_PLACEHOLDER \
    --model_dir AF3_WEIGHT_DIR_PLACEHOLDER \
    --db_dir AF3_DB_DIR_PLACEHOLDER \
    --jackhmmer_binary_path HMMER_DIR_PLACEHOLDER/jackhmmer \
    --hmmalign_binary_path HMMER_DIR_PLACEHOLDER/hmmalign \
    --hmmbuild_binary_path HMMER_DIR_PLACEHOLDER/hmmbuild \
    --hmmsearch_binary_path HMMER_DIR_PLACEHOLDER/hmmsearch \
    --nhmmer_binary_path HMMER_DIR_PLACEHOLDER/nhmmer \
    --output_dir PREFOLD_DIR_PLACEHOLDER

if [[ -f "PREFOLD_DIR_PLACEHOLDER/PROTEIN_LOWER_PLACEHOLDER/PROTEIN_LOWER_PLACEHOLDER_summary_confidences.json" ]]; then
    touch "TOKEN_PLACEHOLDER"
else
    echo "Expected AF3 output not found: PREFOLD_DIR_PLACEHOLDER/PROTEIN_LOWER_PLACEHOLDER/PROTEIN_LOWER_PLACEHOLDER_summary_confidences.json" >&2
    exit 1
fi

echo "Job completed at: $(date)"
EOFSCRIPT

    # Replace placeholders
    sed -i "s|PROTEIN_LOWER_PLACEHOLDER|${protein_lower}|g" "$JOB_SCRIPT"
    sed -i "s|PROTEIN_PLACEHOLDER|${protein}|g" "$JOB_SCRIPT"
    sed -i "s|CONSTRAINT_PLACEHOLDER|${CONSTRAINT}|g" "$JOB_SCRIPT"
    sed -i "s|TIME_LIMIT_PLACEHOLDER|${TIME_LIMIT}|g" "$JOB_SCRIPT"
    sed -i "s|MEMORY_PLACEHOLDER|${MEMORY}|g" "$JOB_SCRIPT"
    sed -i "s|CPUS_PLACEHOLDER|${CPUS_PER_TASK}|g" "$JOB_SCRIPT"
    sed -i "s|LOG_DIR_PLACEHOLDER|${LOG_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|CONDA_ENV_PLACEHOLDER|${CONDA_ENV}|g" "$JOB_SCRIPT"
    sed -i "s|AF3_EXE_PLACEHOLDER|${AF3_EXE}|g" "$JOB_SCRIPT"
    sed -i "s|INPUT_JSON_PLACEHOLDER|${INPUT_JSON}|g" "$JOB_SCRIPT"
    sed -i "s|AF3_WEIGHT_DIR_PLACEHOLDER|${AF3_WEIGHT_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|AF3_DB_DIR_PLACEHOLDER|${AF3_DB_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|HMMER_DIR_PLACEHOLDER|${HMMER_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|PREFOLD_DIR_PLACEHOLDER|${PREFOLD_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|TOKEN_PLACEHOLDER|${TOKEN}|g" "$JOB_SCRIPT"

    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted AF3 prefold for ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${PREFOLD_DIR}/submitted_jobs.txt"
        return 0
    else
        log_error "Failed to submit AF3 prefold for ${protein}"
        return 1
    fi
}

# ============================================================================
# Main
# ============================================================================

log_info "Starting AF3 prefold job submission"
log_info "Task root: ${TASK_ROOT}"

PROTEINS=$(get_proteins)
log_info "Processing proteins: ${PROTEINS}"

SUBMITTED=0

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: ${PROTEIN}"

    if submit_prefold_job "$PROTEIN"; then
        SUBMITTED=$((SUBMITTED + 1))
    fi
done

log_info "Submitted ${SUBMITTED} AF3 prefold jobs"
log_info "Monitor jobs with: squeue -u \$USER"
