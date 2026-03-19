#!/bin/bash
#
# AF3 Batch Prediction Automation Script
# This script submits multiple sbatch jobs for AlphaFold3 predictions
# Based on the Snakemake workflow, focusing only on the prediction step
#

set -e

# ============================================================================
# Configuration - Update these paths as needed
# ============================================================================

# Source conda configuration
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

# Load configuration from environment or use defaults
TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
SCRIPT_ROOT="/home/yangl_pacagen_com/screening_workflow/scripts"

# Tools
AF3_EXE="/home/yangl_pacagen_com/Applications/alphafold3/run_alphafold.py"

# Weight and database directories
AF3_WEIGHT_DIR="${MASTER_AF3_WEIGHT_DIR:-/home/yangl_pacagen_com/Applications/model_weights/AF3}"
AF3_DB_DIR="${MASTER_AF3_DB_DIR:-/home/yangl_pacagen_com/Applications/model_weights/af3_db}"

# Number of batches to split jobs into
N_BATCHES=40

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

# Function to get protein list from environment or config
get_proteins() {
    if [ -n "$MASTER_PROTEINS" ]; then
        echo "$MASTER_PROTEINS"
    else
        # Default protein list if not set by master script
        echo "JAK1JH1"
    fi
}

# Function to submit a single AF3 batch job
submit_af3_batch() {
    local protein=$1
    local batch_id=$2

    local INPUT_DIR="${TASK_ROOT}/${protein}/fine_screening/AF3/input"
    local LOG_DIR="${TASK_ROOT}/${protein}/fine_screening/AF3/logs"
    local FINE_DIR="${TASK_ROOT}/${protein}/fine_screening"
    local OUTPUT_DIR="${FINE_DIR}/AF3/output"
    local TOKEN_DIR="${OUTPUT_DIR}/token"
    local SELECTED_CSV="${TASK_ROOT}/${protein}/initial_screening/selected.csv"
    local LEDGER_FILE="${OUTPUT_DIR}/job_results.tsv"

    # Check if input directory exists
    if [ ! -d "$INPUT_DIR" ]; then
        log_error "Input directory not found: $INPUT_DIR"
        return 1
    fi

    # Check if selected CSV exists
    if [ ! -f "$SELECTED_CSV" ]; then
        log_error "Selected CSV not found: $SELECTED_CSV"
        return 1
    fi

    # Create output directories
    mkdir -p "$TOKEN_DIR"
    mkdir -p "$LOG_DIR"

    # Create SLURM job script
    local JOB_SCRIPT="${LOG_DIR}/slurm_batch_${batch_id}.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=af3_${protein}_b${batch_id}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --mem=${MEMORY}
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --partition=${PARTITION}
#SBATCH --output=${LOG_DIR}/slurm_batch_${batch_id}_%j.out
#SBATCH --error=${LOG_DIR}/slurm_batch_${batch_id}_%j.err

# Error handling
set -e

# Log start time
echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Processing batch ${batch_id} for protein ${protein}"

# Activate conda environment
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate af3_test
source "${SCRIPT_ROOT}/pipeline/lib/job_result_ledger.sh"

# Set AF3-specific GPU environment variables
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

# Calculate job range for this batch
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

# Create batch-level temporary directory on local NVMe/SSD
BATCH_LOCAL_OUT=/tmp/af3_batch_${protein}_${batch_id}_\${SLURM_JOB_ID}
mkdir -p \$BATCH_LOCAL_OUT

# Process each job in the batch
SUCCESS_COUNT=0
FAIL_COUNT=0

for i in \$(seq \$START_IDX \$((\$END_IDX - 1))); do
    JSON_FILE="${INPUT_DIR}/${protein}_\${i}.json"

    if [ ! -f "\$JSON_FILE" ]; then
        echo "Warning: JSON file not found: \$JSON_FILE"
        FAIL_COUNT=\$((FAIL_COUNT + 1))
        continue
    fi

    # Use local NVMe for temporary storage - subdirectory per job
    LOCAL_OUT=\$BATCH_LOCAL_OUT/job_\${i}
    mkdir -p \$LOCAL_OUT

    echo "Processing job \$i (batch ${batch_id})"

    # Run AF3 prediction
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

# Copy results directly to output directory
echo "Copying results for batch ${batch_id}..."

for job_dir in \$BATCH_LOCAL_OUT/job_*; do
    if [ ! -d "\$job_dir" ]; then
        continue
    fi

    job_num=\$(basename "\$job_dir" | sed 's/job_//')

    # Find the compound directory (AF3 creates a subdirectory with compound name)
    compound_dir=\$(find "\$job_dir" -mindepth 1 -maxdepth 1 -type d | head -1)

    if [ -z "\$compound_dir" ] || [ ! -d "\$compound_dir" ]; then
        echo "Warning: No compound directory found in \$job_dir"
        continue
    fi

    job_output="${OUTPUT_DIR}/job_\${job_num}"
    mkdir -p "\$job_output"

    # Extract and save SMILES (but don't copy the data.json file)
    data_json=\$(find "\$compound_dir" -name "*_data.json" | head -1)
    if [ -f "\$data_json" ]; then
        python3 -c "
import json
with open('\$data_json') as f:
    data = json.load(f)
smiles = data['sequences'][-1]['ligand']['smiles']
with open('\$job_output/smiles.txt', 'w') as f:
    f.write(smiles)
" 2>/dev/null || echo "Warning: Could not extract SMILES for job \$job_num"
    fi

    # Copy only the confidence files (actual prediction scores)
    cp "\$compound_dir"/*_summary_confidences.json "\$job_output/" 2>/dev/null || true
    cp "\$compound_dir"/*_confidences.json "\$job_output/" 2>/dev/null || true
    cp "\$compound_dir"/*_model.cif "\$job_output/" 2>/dev/null || true

    echo "  Copied results for job \$job_num"
done

# Create per-compound completion tokens (for Snakemake compatibility)
for i in \$(seq \$START_IDX \$((\$END_IDX - 1))); do
    touch "${TOKEN_DIR}/\${i}.done"
done

# Validate copied outputs and append batch result to job ledger
EXPECTED_COUNT=\$((END_IDX - START_IDX))
COMPLETED_COUNT=0
for i in \$(seq \$START_IDX \$((\$END_IDX - 1))); do
    if compgen -G "${OUTPUT_DIR}/job_\${i}/*_summary_confidences.json" > /dev/null; then
        COMPLETED_COUNT=\$((COMPLETED_COUNT + 1))
    fi
done

JOB_STATUS="failed"
if [ "\$EXPECTED_COUNT" -gt 0 ] && [ "\$COMPLETED_COUNT" -eq "\$EXPECTED_COUNT" ]; then
    JOB_STATUS="success"
fi

append_job_result "${LEDGER_FILE}" "\${SLURM_JOB_ID:-unknown}" "${protein}" "batch_${batch_id}" "\$JOB_STATUS" "\$COMPLETED_COUNT" "\$EXPECTED_COUNT"
echo "Batch ledger status: \${JOB_STATUS} (\${COMPLETED_COUNT}/\${EXPECTED_COUNT})"

# Clean up local temporary directory
rm -rf \$BATCH_LOCAL_OUT

echo "Job completed at: \$(date)"
echo "Results copied to: ${OUTPUT_DIR}/"
EOF

    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted batch ${batch_id} for protein ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${OUTPUT_DIR}/job_ids.txt"
        return 0
    else
        log_error "Failed to submit batch ${batch_id} for protein ${protein}"
        return 1
    fi
}

# ============================================================================
# Main Script
# ============================================================================

log_info "Starting AF3 batch job submission"
log_info "Number of batches per protein: ${N_BATCHES}"

# Get list of proteins to process
PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    # Check if input token exists (prerequisite)
    INPUT_TOKEN="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/input/af3_input.done"
    if [ ! -f "$INPUT_TOKEN" ]; then
        log_error "Input token not found for ${PROTEIN}: ${INPUT_TOKEN}"
        log_error "Please run the input preparation step first"
        continue
    fi

    # Get number of SMILES to process
    SELECTED_CSV="${TASK_ROOT}/${PROTEIN}/initial_screening/selected.csv"
    if [ ! -f "$SELECTED_CSV" ]; then
        log_error "Selected CSV not found: ${SELECTED_CSV}"
        continue
    fi

    N_SMILES=$(tail -n +2 "$SELECTED_CSV" | wc -l)
    BATCH_SIZE=$(( ($N_SMILES + $N_BATCHES - 1) / $N_BATCHES ))

    log_info "Total SMILES for ${PROTEIN}: ${N_SMILES}"
    log_info "Batch size: ${BATCH_SIZE}"

    # Clear previous job list
    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/output"
    mkdir -p "$OUTPUT_DIR"
    > "${OUTPUT_DIR}/job_ids.txt"
    > "${OUTPUT_DIR}/job_results.tsv"

    # Submit jobs for each batch
    SUBMITTED=0
    for BATCH_ID in $(seq 0 $((N_BATCHES - 1))); do
        # Check if this batch has jobs
        START_IDX=$((BATCH_ID * BATCH_SIZE))
        if [ $START_IDX -lt $N_SMILES ]; then
            if submit_af3_batch "$PROTEIN" "$BATCH_ID"; then
                SUBMITTED=$((SUBMITTED + 1))
            fi
        else
            log_info "Batch ${BATCH_ID} has no jobs, skipping"
        fi
    done

    log_info "Submitted ${SUBMITTED} batch jobs for protein ${PROTEIN}"
    log_info "Job IDs saved to: ${OUTPUT_DIR}/submitted_jobs.txt"
done

log_info "All batch jobs submitted!"
log_info ""
log_info "Monitor jobs with: squeue -u \$USER"
log_info "Check job outputs in: \${TASK_ROOT}/\${PROTEIN}/fine_screening/AF3/output/"
