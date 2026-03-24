#!/bin/bash
#
# Boltz2 Batch Prediction Automation Script
# This script submits multiple sbatch jobs for Boltz2 predictions
# Based on the Snakemake workflow, focusing only on the prediction step
#

set -e

# ============================================================================
# Configuration - Update these paths as needed
# ============================================================================

# Source conda configuration
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

# Load configuration from environment or use defaults
TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Tools
BOLTZ2_EXE="/home/ubuntu/miniconda3/envs/boltz/bin/boltz"

# Number of batches to split jobs into
N_BATCHES=40

# SLURM configuration
# PARTITION="gpu"           # Change to your GPU partition name
TIME_LIMIT="48:00:00"     # 4 hours per job
MEMORY="15G"              # Memory per job
CPUS_PER_TASK=4
# GPUS_PER_NODE=1

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

get_selected_csv() {
    local protein=$1
    local rel_path="${MASTER_SELECTED_REL_PATH:-initial_screening/selected.csv}"
    echo "${TASK_ROOT}/${protein}/${rel_path}"
}

# Function to submit a single Boltz2 batch job
submit_boltz2_batch() {
    local protein=$1
    local batch_id=$2

    local INPUT_DIR="/home/ubuntu/${protein}/boltz2_tmp/input"
    local LOG_DIR="/home/ubuntu/${protein}/boltz2_tmp/logs"
    local FINE_DIR="${TASK_ROOT}/${protein}/fine_screening"
    local OUTPUT_DIR="${FINE_DIR}/Boltz2/output"
    local TOKEN_DIR="${OUTPUT_DIR}/token"
    local SELECTED_CSV="$(get_selected_csv "${protein}")"

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
#SBATCH --job-name=boltz2_${protein}_b${batch_id}
#SBATCH --constraint=g5.xlarge
#SBATCH --time=${TIME_LIMIT}
#SBATCH --mem=${MEMORY}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
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
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate boltz

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
BATCH_LOCAL_OUT=/tmp/boltz2_batch_${protein}_${batch_id}_\${SLURM_JOB_ID}
mkdir -p \$BATCH_LOCAL_OUT

# Process each job in the batch
SUCCESS_COUNT=0
FAIL_COUNT=0

for i in \$(seq \$START_IDX \$((\$END_IDX - 1))); do
    YAML_FILE="${INPUT_DIR}/\${i}.yaml"

    if [ ! -f "\$YAML_FILE" ]; then
        echo "Warning: YAML file not found: \$YAML_FILE"
        FAIL_COUNT=\$((FAIL_COUNT + 1))
        continue
    fi

    # Use local NVMe for temporary storage - subdirectory per job
    LOCAL_OUT=\$BATCH_LOCAL_OUT/job_\${i}
    mkdir -p \$LOCAL_OUT

    echo "Processing job \$i (batch ${batch_id})"

    # Run Boltz2 prediction
    if "${BOLTZ2_EXE}" predict "\$YAML_FILE" \\
        --out_dir=\$LOCAL_OUT \\
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

# Compress results
echo "Compressing results for batch ${batch_id}..."
cd \$BATCH_LOCAL_OUT

# Find and compress only the relevant output files (affinity, confidence, models)
find . \\( -name "affinity_*.json" -o -name "confidence*.json" -o -name "*model*.cif" \\) -print0 | \\
    tar -czf batch_${batch_id}.tar.gz --null -T -

# Copy tar file to final output directory
cp batch_${batch_id}.tar.gz "${OUTPUT_DIR}/"

# Create completion token
touch "${TOKEN_DIR}/batch_${batch_id}.done"

# Clean up local temporary directory
rm -rf \$BATCH_LOCAL_OUT

echo "Job completed at: \$(date)"
echo "Results saved to: ${OUTPUT_DIR}/batch_${batch_id}.tar.gz"
EOF

    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted batch ${batch_id} for protein ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${OUTPUT_DIR}/submitted_jobs.txt"
        return 0
    else
        log_error "Failed to submit batch ${batch_id} for protein ${protein}"
        return 1
    fi
}

# ============================================================================
# Main Script
# ============================================================================

log_info "Starting Boltz2 batch job submission"
log_info "Number of batches per protein: ${N_BATCHES}"

# Get list of proteins to process
PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    # Check if input token exists (prerequisite)
    INPUT_TOKEN="/home/ubuntu/${PROTEIN}/boltz2_tmp/boltz_input.done"
    if [ ! -f "$INPUT_TOKEN" ]; then
        log_error "Input token not found for ${PROTEIN}: ${INPUT_TOKEN}"
        log_error "Please run the input preparation step first"
        continue
    fi

    # Get number of SMILES to process
    SELECTED_CSV="$(get_selected_csv "${PROTEIN}")"
    if [ ! -f "$SELECTED_CSV" ]; then
        log_error "Selected CSV not found: ${SELECTED_CSV}"
        continue
    fi

    N_SMILES=$(tail -n +2 "$SELECTED_CSV" | wc -l)
    BATCH_SIZE=$(( ($N_SMILES + $N_BATCHES - 1) / $N_BATCHES ))

    log_info "Total SMILES for ${PROTEIN}: ${N_SMILES}"
    log_info "Batch size: ${BATCH_SIZE}"

    # Clear previous job list
    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/output"
    mkdir -p "$OUTPUT_DIR"
    > "${OUTPUT_DIR}/submitted_jobs.txt"

    # Submit jobs for each batch
    SUBMITTED=0
    for BATCH_ID in $(seq 0 $((N_BATCHES - 1))); do
        # Check if this batch has jobs
        START_IDX=$((BATCH_ID * BATCH_SIZE))
        if [ $START_IDX -lt $N_SMILES ]; then
            if submit_boltz2_batch "$PROTEIN" "$BATCH_ID"; then
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
log_info "Check job outputs in: \${TASK_ROOT}/\${PROTEIN}/fine_screening/Boltz2/output/"
