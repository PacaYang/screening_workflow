#!/bin/bash
#
# Standalone script for Snakemake checkpoint: write_af3_input
# Generates AF3 JSON inputs with ligands for fine screening.
# Submits one SLURM job per protein and waits for all to complete.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
EXE="/home/yangl_pacagen_com/screening_workflow/scripts/input/gen_af3_json_with_cmpds.py"

PARTITION="g24"
MEMORY="8G"
CPUS=2
TIME_LIMIT="2:00:00"
CONDA_ENV="HMSA_test"
POLL_INTERVAL=30

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

# ============================================================================
# Main
# ============================================================================

log_info "Starting write_af3_input (SLURM mode)"

PROTEINS=$(get_proteins)

declare -A JOB_IDS   # protein -> job_id
declare -A TOKENS    # protein -> token path

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    OUTDIR="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/input"
    TOKEN="${OUTDIR}/af3_input.done"
    PREFOLD_JSON="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/prefold/${PROTEIN}/${PROTEIN}_data.json"
    SELECTED="${TASK_ROOT}/${PROTEIN}/initial_screening/selected.csv"
    PREFOLD_TOKEN="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/prefold/prefold.done"
    LOG_DIR="${OUTDIR}/logs"

    TOKENS[$PROTEIN]="$TOKEN"

    if [ -f "$TOKEN" ]; then
        log_info "Already completed for ${PROTEIN}, skipping"
        continue
    fi

    if [ ! -f "$PREFOLD_TOKEN" ]; then
        log_error "prefold.done not found for ${PROTEIN}: $PREFOLD_TOKEN"
        exit 1
    fi

    if [ ! -f "$SELECTED" ]; then
        log_error "selected.csv not found for ${PROTEIN}: $SELECTED"
        exit 1
    fi

    mkdir -p "$OUTDIR" "$LOG_DIR"

    JOB_SCRIPT="${LOG_DIR}/slurm_write_af3_input.sh"

    cat > "$JOB_SCRIPT" <<SLURM_SCRIPT
#!/bin/bash
#SBATCH --job-name=write_af3_${PROTEIN}
#SBATCH --partition=${PARTITION}
#SBATCH --mem=${MEMORY}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOG_DIR}/write_af3_input_%j.out
#SBATCH --error=${LOG_DIR}/write_af3_input_%j.err

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

echo "[\\$(date '+%Y-%m-%d %H:%M:%S')] Writing AF3 inputs for ${PROTEIN}"
python ${EXE} \\
    --output-dir ${OUTDIR} \\
    --input-json ${PREFOLD_JSON} \\
    --smiles-file ${SELECTED} \\
    --smiles-col SMILES

touch ${TOKEN}
echo "[\\$(date '+%Y-%m-%d %H:%M:%S')] Done: ${TOKEN}"
SLURM_SCRIPT

    chmod +x "$JOB_SCRIPT"

    if ! JOB_ID=$(sbatch --parsable "$JOB_SCRIPT"); then
        log_error "sbatch submission failed for ${PROTEIN}"
        exit 1
    fi
    log_info "Submitted SLURM job for ${PROTEIN}: job ID ${JOB_ID}"
    JOB_IDS[$PROTEIN]="$JOB_ID"
done

# ============================================================================
# Wait for all submitted jobs to finish
# ============================================================================

if [ ${#JOB_IDS[@]} -eq 0 ]; then
    log_info "No jobs submitted (all proteins already done)"
else
    log_info "Waiting for ${#JOB_IDS[@]} SLURM job(s) to complete..."

    while true; do
        ALL_DONE=true
        for PROTEIN in "${!JOB_IDS[@]}"; do
            JOB_ID="${JOB_IDS[$PROTEIN]}"
            if squeue -j "$JOB_ID" -h >/dev/null 2>&1 && [ -n "$(squeue -j "$JOB_ID" -h 2>/dev/null)" ]; then
                ALL_DONE=false
                log_info "Job ${JOB_ID} (${PROTEIN}) still running..."
            fi
        done
        if $ALL_DONE; then
            break
        fi
        sleep "$POLL_INTERVAL"
    done

    log_info "All jobs finished. Verifying outputs..."

    FAILED=0
    for PROTEIN in "${!JOB_IDS[@]}"; do
        TOKEN="${TOKENS[$PROTEIN]}"
        if [ ! -f "$TOKEN" ]; then
            log_error "af3_input.done not found for ${PROTEIN} after job completion: $TOKEN"
            FAILED=1
        else
            log_info "Verified: ${PROTEIN} af3_input.done exists"
        fi
    done

    if [ "$FAILED" -eq 1 ]; then
        log_error "One or more proteins failed to produce af3_input.done"
        exit 1
    fi
fi

log_info "Done"
