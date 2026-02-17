#!/bin/bash
#
# Vina Batch Docking Automation Script
# This script submits multiple sbatch jobs for Vina docking
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
SCRIPT_ROOT="/home/ubuntu/screening_workflow/scripts"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"

# Tools
VINA_EXE="/home/ubuntu/screening_workflow/scripts/docking.py"

# SLURM configuration
TIME_LIMIT="48:00:00"     # 8 hours per job
MEMORY="15G"              # Memory per job
CPUS_PER_TASK=2

# GPU configuration (Vina might not need GPU, adjust as needed)
# Set to empty string if no GPU needed
GPU_REQUEST=""            # Empty = no GPU request
# GPU_REQUEST="--gres=gpu:1"  # Uncomment if GPU is needed

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

# Function to extract docking box from sequences.csv
get_box_params() {
    local protein=$1
    local seq_file=$2

    # Use Python to parse the docking box and output as JSON
    python3 <<EOF
import pandas as pd
import ast
import json

df = pd.read_csv("${seq_file}")
cell = df.loc[df['name'] == "${protein}", 'docking box'].iloc[0]

# Parse the box specification
if isinstance(cell, str):
    vals = ast.literal_eval(cell)
else:
    vals = list(cell)

# Normalize flat list to nested: [6] -> [[6]]
if vals and not isinstance(vals[0], (list, tuple)):
    vals = [vals]

print(json.dumps(vals))
EOF
}

# Function to submit a single Vina docking job
submit_vina_job() {
    local protein=$1
    local part=$2

    local FINE_DIR="${TASK_ROOT}/${protein}/fine_screening"
    local INPUT_CSV="${FINE_DIR}/Vina/input/${part}.csv"
    local OUTPUT_DIR="${FINE_DIR}/Vina/output/${part}"
    local TOKEN_FILE="${FINE_DIR}/Vina/output/${part}.done"
    local PDB_FILE="${TASK_ROOT}/Input/protein_file/${protein}/${protein}.pdb"

    # Check if input CSV exists
    if [ ! -f "$INPUT_CSV" ]; then
        log_error "Input CSV not found: $INPUT_CSV"
        return 1
    fi

    # Check if PDB file exists
    if [ ! -f "$PDB_FILE" ]; then
        log_error "PDB file not found: $PDB_FILE"
        return 1
    fi

    # Get docking box parameters as JSON
    BOXES_JSON=$(get_box_params "$protein" "$SEQS_CSV")

    if [ -z "$BOXES_JSON" ]; then
        log_error "Failed to extract docking box parameters for ${protein}"
        return 1
    fi

    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$(dirname "$TOKEN_FILE")"

    # Create SLURM job script
    local JOB_SCRIPT="${FINE_DIR}/Vina/output/slurm_${part}.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=vina_${protein}_${part}
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
#SBATCH --output=${FINE_DIR}/Vina/output/slurm_${part}_%j.out
#SBATCH --error=${FINE_DIR}/Vina/output/slurm_${part}_%j.err

# Error handling
set -e

# Log start time
echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Processing part ${part} for protein ${protein}"
echo "Docking boxes: ${BOXES_JSON}"

# Activate conda environment
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate vina_new

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Run Vina docking
echo "Starting Vina docking..."
python "${VINA_EXE}" \\
    --smiles "${INPUT_CSV}" \\
    --pdb "${PDB_FILE}" \\
    --boxes '${BOXES_JSON}' \\
    --output "${OUTPUT_DIR}" \\
    --smiles-col "ligand_description" \\
    --skip-docked

# Check if docking completed successfully
if [ \$? -eq 0 ]; then
    echo "Vina docking completed successfully"
    touch "${TOKEN_FILE}"
else
    echo "Vina docking failed"
    exit 1
fi

echo "Job completed at: \$(date)"
EOF

    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted part ${part} for protein ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${FINE_DIR}/Vina/output/submitted_jobs.txt"
        return 0
    else
        log_error "Failed to submit part ${part} for protein ${protein}"
        return 1
    fi
}

# ============================================================================
# Main Script
# ============================================================================

log_info "Starting Vina batch job submission"

# Get list of proteins to process
PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    # Check if input directory exists
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
    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/output"
    mkdir -p "$OUTPUT_DIR"
    > "${OUTPUT_DIR}/submitted_jobs.txt"

    # Submit jobs for each part
    SUBMITTED=0
    for INPUT_FILE in "${INPUT_FILES[@]}"; do
        # Extract part name (e.g., "part_0" from "/path/to/part_0.csv")
        PART=$(basename "$INPUT_FILE" .csv)

        if submit_vina_job "$PROTEIN" "$PART"; then
            SUBMITTED=$((SUBMITTED + 1))
        fi
    done

    log_info "Submitted ${SUBMITTED} docking jobs for protein ${PROTEIN}"
    log_info "Job IDs saved to: ${OUTPUT_DIR}/submitted_jobs.txt"
done

log_info "All docking jobs submitted!"
log_info ""
log_info "Monitor jobs with: squeue -u \$USER"
log_info "Check job outputs in: \${TASK_ROOT}/\${PROTEIN}/fine_screening/Vina/output/"
