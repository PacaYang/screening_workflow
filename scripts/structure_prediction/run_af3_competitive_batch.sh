#!/bin/bash
#
# AF3 Competitive Folding Batch Prediction Script
# Submits SLURM jobs for AlphaFold3 predictions with 2 competing ligands (Z and Y).
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

AF3_EXE="/home/ubuntu/Applications/alphafold3/run_alphafold.py"
AF3_WEIGHT_DIR="/shared/programs/af3_weights"
AF3_DB_DIR="/shared/programs/af3_data"

N_BATCHES=40

TIME_LIMIT="48:00:00"
MEMORY="15G"
CPUS_PER_TASK=4

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

get_selected_csv() {
    local protein=$1
    local rel_path="${MASTER_SELECTED_REL_PATH:-initial_screening/selected.csv}"
    echo "${TASK_ROOT}/${protein}/${rel_path}"
}

submit_competitive_batch() {
    local protein=$1
    local batch_id=$2

    local PROTEIN_LOWER=$(echo "$protein" | tr '[:upper:]' '[:lower:]')
    local INPUT_DIR="${TASK_ROOT}/${protein}/fine_screening/AF3_competitive/input"
    local LOG_DIR="${TASK_ROOT}/${protein}/fine_screening/AF3_competitive/logs"
    local OUTPUT_DIR="${TASK_ROOT}/${protein}/fine_screening/AF3_competitive/output"
    local TOKEN_DIR="${OUTPUT_DIR}/token"
    local SELECTED_CSV="$(get_selected_csv "${protein}")"

    if [ ! -d "$INPUT_DIR" ]; then
        log_error "Input directory not found: $INPUT_DIR"
        return 1
    fi

    if [ ! -f "$SELECTED_CSV" ]; then
        log_error "Selected CSV not found: $SELECTED_CSV"
        return 1
    fi

    mkdir -p "$TOKEN_DIR"
    mkdir -p "$LOG_DIR"

    local JOB_SCRIPT="${LOG_DIR}/slurm_batch_${batch_id}.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=af3_comp_${protein}_b${batch_id}
#SBATCH --constraint=g5.xlarge
#SBATCH --time=${TIME_LIMIT}
#SBATCH --mem=${MEMORY}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --output=${LOG_DIR}/slurm_batch_${batch_id}_%j.out
#SBATCH --error=${LOG_DIR}/slurm_batch_${batch_id}_%j.err

set -e

echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Competitive folding batch ${batch_id} for protein ${protein}"

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate af3

export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

SMILES_FILE="${SELECTED_CSV}"
N_SMILES=\$(tail -n +2 "\$SMILES_FILE" | wc -l)

BATCH_SIZE=\$(( (\$N_SMILES + ${N_BATCHES} - 1) / ${N_BATCHES} ))
START_IDX=\$(( ${batch_id} * \$BATCH_SIZE ))
END_IDX=\$(( \$START_IDX + \$BATCH_SIZE ))
if [ \$END_IDX -gt \$N_SMILES ]; then
    END_IDX=\$N_SMILES
fi

echo "Total SMILES: \$N_SMILES"
echo "Batch size: \$BATCH_SIZE"
echo "Processing jobs \$START_IDX to \$((\$END_IDX - 1))"

BATCH_LOCAL_OUT=/tmp/af3_comp_${protein}_${batch_id}_\${SLURM_JOB_ID}
mkdir -p \$BATCH_LOCAL_OUT

SUCCESS_COUNT=0
FAIL_COUNT=0

PROTEIN_LOWER="${PROTEIN_LOWER}"

for i in \$(seq \$START_IDX \$((\$END_IDX - 1))); do
    JSON_FILE="${INPUT_DIR}/\${PROTEIN_LOWER}_\${i}.json"

    if [ ! -f "\$JSON_FILE" ]; then
        echo "Warning: JSON file not found: \$JSON_FILE"
        FAIL_COUNT=\$((FAIL_COUNT + 1))
        continue
    fi

    LOCAL_OUT=\$BATCH_LOCAL_OUT/job_\${i}
    mkdir -p \$LOCAL_OUT

    echo "Processing competitive job \$i (batch ${batch_id})"

    if python ${AF3_EXE} \\
        --json_path=\$JSON_FILE \\
        --model_dir=${AF3_WEIGHT_DIR} \\
        --db_dir=${AF3_DB_DIR} \\
        --output_dir=\$LOCAL_OUT \\
        --norun_data_pipeline \\
        ; then
        SUCCESS_COUNT=\$((SUCCESS_COUNT + 1))
        echo "  ✓ Job \$i completed successfully"
    else
        FAIL_COUNT=\$((FAIL_COUNT + 1))
        echo "  ✗ Job \$i failed"
    fi
done

echo ""
echo "Batch ${batch_id} summary:"
echo "  Successful: \$SUCCESS_COUNT"
echo "  Failed: \$FAIL_COUNT"
echo ""

echo "Compressing results for batch ${batch_id}..."
cd \$BATCH_LOCAL_OUT

find . \\( -name "*_data.json" -o -name "*_summary_confidences.json" -o -name "*_confidences.json" -o -name "*_model.cif" \\) -print0 | \\
    tar -czf batch_${batch_id}.tar.gz --null -T -

cp batch_${batch_id}.tar.gz "${OUTPUT_DIR}/"

for i in \$(seq \$START_IDX \$((\$END_IDX - 1))); do
    touch "${TOKEN_DIR}/\${i}.done"
done

rm -rf \$BATCH_LOCAL_OUT

echo "Job completed at: \$(date)"
echo "Results saved to: ${OUTPUT_DIR}/batch_${batch_id}.tar.gz"
EOF

    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted competitive batch ${batch_id} for protein ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${OUTPUT_DIR}/submitted_jobs.txt"
        return 0
    else
        log_error "Failed to submit competitive batch ${batch_id} for protein ${protein}"
        return 1
    fi
}

# ============================================================================
# Main
# ============================================================================

log_info "Starting AF3 competitive folding batch submission"
log_info "Number of batches per protein: ${N_BATCHES}"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    INPUT_TOKEN="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3_competitive/input/af3_competitive_input.done"
    if [ ! -f "$INPUT_TOKEN" ]; then
        log_error "Input token not found for ${PROTEIN}: ${INPUT_TOKEN}"
        log_error "Please run run_write_af3_competitive_input.sh first"
        continue
    fi

    SELECTED_CSV="$(get_selected_csv "${PROTEIN}")"
    if [ ! -f "$SELECTED_CSV" ]; then
        log_error "Selected CSV not found: ${SELECTED_CSV}"
        continue
    fi

    N_SMILES=$(tail -n +2 "$SELECTED_CSV" | wc -l)
    BATCH_SIZE=$(( ($N_SMILES + $N_BATCHES - 1) / $N_BATCHES ))

    log_info "Total SMILES for ${PROTEIN}: ${N_SMILES}"
    log_info "Batch size: ${BATCH_SIZE}"

    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3_competitive/output"
    mkdir -p "$OUTPUT_DIR"
    > "${OUTPUT_DIR}/submitted_jobs.txt"

    SUBMITTED=0
    for BATCH_ID in $(seq 0 $((N_BATCHES - 1))); do
        START_IDX=$((BATCH_ID * BATCH_SIZE))
        if [ $START_IDX -lt $N_SMILES ]; then
            if submit_competitive_batch "$PROTEIN" "$BATCH_ID"; then
                SUBMITTED=$((SUBMITTED + 1))
            fi
        else
            log_info "Batch ${BATCH_ID} has no jobs, skipping"
        fi
    done

    log_info "Submitted ${SUBMITTED} competitive batch jobs for protein ${PROTEIN}"
    log_info "Job IDs saved to: ${OUTPUT_DIR}/submitted_jobs.txt"
done

log_info "All competitive folding batch jobs submitted!"
log_info "Monitor jobs with: squeue -u \$USER"
