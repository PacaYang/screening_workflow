#!/bin/bash
#
# Standalone script for Snakemake rule: prefold_boltz2
# Submits a SLURM job per protein to run Boltz2 protein-only folding.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
BOLTZ_EXE="/home/ubuntu/miniconda3/envs/boltz/bin/boltz"
CONDA_ENV="boltz"

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

    local PREFOLD_DIR="${TASK_ROOT}/${protein}/fine_screening/Boltz2/prefold"
    local INPUT_YAML="${PREFOLD_DIR}/${protein}.yaml"
    local OUTDIR="/home/ubuntu/${protein}/boltz2_tmp"
    local CONFIDENCE="${OUTDIR}/boltz_results_${protein}/predictions/${protein}/confidence_${protein}_model_0.json"
    local LOG_DIR="${PREFOLD_DIR}/logs"

    if [ -f "$CONFIDENCE" ]; then
        log_info "Confidence JSON already exists for ${protein}, skipping"
        return 0
    fi

    if [ ! -f "$INPUT_YAML" ]; then
        log_error "Input YAML not found for ${protein}: ${INPUT_YAML}"
        log_error "Please run run_write_prefold_boltz2.sh first"
        return 1
    fi

    mkdir -p "$LOG_DIR"

    local JOB_SCRIPT="${LOG_DIR}/slurm_prefold.sh"

    cat > "$JOB_SCRIPT" <<'EOFSCRIPT'
#!/bin/bash
#SBATCH --job-name=boltz2_prefold_PROTEIN_PLACEHOLDER
#SBATCH --constraint=CONSTRAINT_PLACEHOLDER
#SBATCH --time=TIME_LIMIT_PLACEHOLDER
#SBATCH --mem=MEMORY_PLACEHOLDER
#SBATCH --cpus-per-task=CPUS_PLACEHOLDER
#SBATCH --output=LOG_DIR_PLACEHOLDER/slurm_prefold_%j.out
#SBATCH --error=LOG_DIR_PLACEHOLDER/slurm_prefold_%j.err

set -e

echo "Job started at: $(date)"
echo "Running on host: $(hostname)"
echo "Folding protein: PROTEIN_PLACEHOLDER"

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate CONDA_ENV_PLACEHOLDER

BOLTZ_EXE_PLACEHOLDER predict "INPUT_YAML_PLACEHOLDER" \
    --out_dir=OUTDIR_PLACEHOLDER \
    --use_msa_server \
    --override

echo "Job completed at: $(date)"
EOFSCRIPT

    # Replace placeholders
    sed -i "s|PROTEIN_PLACEHOLDER|${protein}|g" "$JOB_SCRIPT"
    sed -i "s|CONSTRAINT_PLACEHOLDER|${CONSTRAINT}|g" "$JOB_SCRIPT"
    sed -i "s|TIME_LIMIT_PLACEHOLDER|${TIME_LIMIT}|g" "$JOB_SCRIPT"
    sed -i "s|MEMORY_PLACEHOLDER|${MEMORY}|g" "$JOB_SCRIPT"
    sed -i "s|CPUS_PLACEHOLDER|${CPUS_PER_TASK}|g" "$JOB_SCRIPT"
    sed -i "s|LOG_DIR_PLACEHOLDER|${LOG_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|CONDA_ENV_PLACEHOLDER|${CONDA_ENV}|g" "$JOB_SCRIPT"
    sed -i "s|BOLTZ_EXE_PLACEHOLDER|${BOLTZ_EXE}|g" "$JOB_SCRIPT"
    sed -i "s|INPUT_YAML_PLACEHOLDER|${INPUT_YAML}|g" "$JOB_SCRIPT"
    sed -i "s|OUTDIR_PLACEHOLDER|${OUTDIR}/|g" "$JOB_SCRIPT"

    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted Boltz2 prefold for ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${PREFOLD_DIR}/submitted_jobs.txt"
        return 0
    else
        log_error "Failed to submit Boltz2 prefold for ${protein}"
        return 1
    fi
}

# ============================================================================
# Main
# ============================================================================

log_info "Starting Boltz2 prefold job submission"
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

log_info "Submitted ${SUBMITTED} Boltz2 prefold jobs"
log_info "Monitor jobs with: squeue -u \$USER"
