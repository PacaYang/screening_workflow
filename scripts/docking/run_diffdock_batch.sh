#!/bin/bash
#
# DiffDock Batch Docking Automation Script
# This script submits multiple sbatch jobs for DiffDock docking
# Based on the Snakemake workflow, focusing only on the docking step
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
DIFFDOCK_DIR="/home/ubuntu/Applications/DiffDock/"
DIFFDOCK_CONFIG="/home/ubuntu/Applications/DiffDock/default_inference_args.yaml"

# SLURM configuration
TIME_LIMIT="48:00:00"     # 12 hours per job (DiffDock can be slow)
MEMORY="15G"              # Memory per job (DiffDock needs more memory)
CPUS_PER_TASK=4

# EC2 instance constraint (if using AWS ParallelCluster)
CONSTRAINT="g5.xlarge"    # Set to empty string if not using constraints
# CONSTRAINT=""

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

# Function to submit a single DiffDock job
submit_diffdock_job() {
    local protein=$1
    local part=$2

    local FINE_DIR="${TASK_ROOT}/${protein}/fine_screening"
    local INPUT_CSV="${FINE_DIR}/Vina/input/${part}.csv"
    local OUTPUT_DIR="${FINE_DIR}/PBSA/DiffDock/output"
    local TOKEN_FILE="${OUTPUT_DIR}/${part}.done"

    # Check if input CSV exists
    if [ ! -f "$INPUT_CSV" ]; then
        log_error "Input CSV not found: $INPUT_CSV"
        return 1
    fi

    # Check if DiffDock directory exists
    if [ ! -d "$DIFFDOCK_DIR" ]; then
        log_error "DiffDock directory not found: $DIFFDOCK_DIR"
        return 1
    fi

    # Check if config file exists
    if [ ! -f "$DIFFDOCK_CONFIG" ]; then
        log_error "DiffDock config not found: $DIFFDOCK_CONFIG"
        return 1
    fi

    # Create output directory
    mkdir -p "$OUTPUT_DIR"

    # Create SLURM job script
    local JOB_SCRIPT="${OUTPUT_DIR}/slurm_${part}.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=diffdock_${protein}_${part}
EOF

    # Add constraint if specified
    if [ -n "$CONSTRAINT" ]; then
        echo "#SBATCH --constraint=${CONSTRAINT}" >> "$JOB_SCRIPT"
    fi

    # Add GPU request if specified
    if [ -n "$GPU_REQUEST" ]; then
        echo "#SBATCH ${GPU_REQUEST}" >> "$JOB_SCRIPT"
    fi

    cat >> "$JOB_SCRIPT" <<EOF
#SBATCH --time=${TIME_LIMIT}
#SBATCH --mem=${MEMORY}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --output=${OUTPUT_DIR}/slurm_${part}_%j.out
#SBATCH --error=${OUTPUT_DIR}/slurm_${part}_%j.err

# Error handling (DiffDock may have some warnings, so using set +eu)
set +eu

# Log start time
echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Processing part ${part} for protein ${protein}"

# Activate conda environment
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate diffdock

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Change to DiffDock directory (required for imports)
cd "${DIFFDOCK_DIR}"

# Run DiffDock inference
echo "Starting DiffDock inference..."
echo "Input CSV: ${INPUT_CSV}"
echo "Config: ${DIFFDOCK_CONFIG}"
echo "Output directory: ${OUTPUT_DIR}"

python -m inference --config "${DIFFDOCK_CONFIG}" \\
    --protein_ligand_csv "${INPUT_CSV}" \\
    --out_dir "${OUTPUT_DIR}"

# Check if DiffDock completed successfully
INFERENCE_EXIT_CODE=\$?

if [ \$INFERENCE_EXIT_CODE -eq 0 ]; then
    echo "DiffDock inference completed successfully"
    touch "${TOKEN_FILE}"
    EXIT_CODE=0
else
    echo "DiffDock inference failed with exit code \$INFERENCE_EXIT_CODE"
    EXIT_CODE=1
fi

echo "Job completed at: \$(date)"
exit \$EXIT_CODE
EOF

    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted part ${part} for protein ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${OUTPUT_DIR}/submitted_jobs.txt"
        return 0
    else
        log_error "Failed to submit part ${part} for protein ${protein}"
        return 1
    fi
}

# ============================================================================
# Main Script
# ============================================================================

log_info "Starting DiffDock batch job submission"

# Get list of proteins to process
PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    # Check if Vina input directory exists (DiffDock uses the same input)
    INPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/input"
    if [ ! -d "$INPUT_DIR" ]; then
        log_error "Input directory not found for ${PROTEIN}: ${INPUT_DIR}"
        log_error "Please run the split_csv preparation step first"
        continue
    fi

    # Find all input CSV files (parts)
    INPUT_FILES=($(ls "$INPUT_DIR"/*.csv 2>/dev/null | sort))

    if [ ${#INPUT_FILES[@]} -eq 0 ]; then
        log_error "No input CSV files found in ${INPUT_DIR}"
        continue
    fi

    N_PARTS=${#INPUT_FILES[@]}
    log_info "Found ${N_PARTS} parts to process for ${PROTEIN}"

    # Clear previous job list
    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/DiffDock/output"
    mkdir -p "$OUTPUT_DIR"
    > "${OUTPUT_DIR}/submitted_jobs.txt"

    # Submit jobs for each part
    SUBMITTED=0
    for INPUT_FILE in "${INPUT_FILES[@]}"; do
        # Extract part name (e.g., "part_0" from "/path/to/part_0.csv")
        PART=$(basename "$INPUT_FILE" .csv)

        if submit_diffdock_job "$PROTEIN" "$PART"; then
            SUBMITTED=$((SUBMITTED + 1))
        fi
    done

    log_info "Submitted ${SUBMITTED} DiffDock jobs for protein ${PROTEIN}"
    log_info "Job IDs saved to: ${OUTPUT_DIR}/submitted_jobs.txt"
done

log_info "All DiffDock jobs submitted!"
log_info ""
log_info "Monitor jobs with: squeue -u \$USER"
log_info "Check job outputs in: \${TASK_ROOT}/\${PROTEIN}/fine_screening/PBSA/DiffDock/output/"
