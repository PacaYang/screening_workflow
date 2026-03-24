#!/bin/bash
#
# Vina Batch Docking Automation Script
# This script prepares receptors, generates split inputs, and submits sbatch jobs.
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
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"

# Tools
VINA_EXE="${SCRIPT_DIR}/docking.py"
SPLIT_CSV_EXE="${SCRIPT_ROOT}/prepare_input/split_csv.py"
PREP_RECEPTOR_SCRIPT="${SCRIPT_ROOT}/prepare_input/prepare_vina_receptor.sh"
VINA_CONDA_ENV="${VINA_CONDA_ENV:-vina_new}"

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

get_selected_csv() {
    local protein=$1
    local rel_path="${MASTER_SELECTED_REL_PATH:-initial_screening/selected.csv}"
    echo "${TASK_ROOT}/${protein}/${rel_path}"
}

# Function to extract docking box from sequences.csv
get_box_params() {
    local protein=$1
    local seq_file=$2

    # Use Python to parse the docking box and output as JSON
    python3 <<EOF_PY
import pandas as pd
import ast
import json


df = pd.read_csv("${seq_file}")
row = df.loc[df['name'] == "${protein}", 'docking box']
if row.empty:
    raise ValueError("Protein not found in sequences.csv: ${protein}")
cell = row.iloc[0]

# Parse the box specification
if isinstance(cell, str):
    vals = ast.literal_eval(cell)
else:
    vals = list(cell)

# Normalize flat list to nested: [6] -> [[6]]
if vals and not isinstance(vals[0], (list, tuple)):
    vals = [vals]

print(json.dumps(vals))
EOF_PY
}

prepare_receptor() {
    local protein=$1
    local receptor_pdbqt="${TASK_ROOT}/${protein}/fine_screening/Vina/receptor/${protein}.pdbqt"

    if ! bash "$PREP_RECEPTOR_SCRIPT" \
        --task-root "$TASK_ROOT" \
        --protein "$protein" \
        --output "$receptor_pdbqt" \
        --conda-env "$VINA_CONDA_ENV" >/dev/null; then
        return 1
    fi

    if [ ! -f "$receptor_pdbqt" ]; then
        log_error "Receptor PDBQT was not created: $receptor_pdbqt"
        return 1
    fi

    echo "$receptor_pdbqt"
    return 0
}

generate_vina_inputs() {
    local protein=$1
    local selected_csv="$(get_selected_csv "${protein}")"
    local pdb_file="${TASK_ROOT}/Input/protein_file/${protein}/${protein}.pdb"
    local input_dir="${TASK_ROOT}/${protein}/fine_screening/Vina/input"

    if [ ! -f "$selected_csv" ]; then
        log_error "Selected compounds file not found for ${protein}: ${selected_csv}"
        return 1
    fi

    if [ ! -f "$pdb_file" ]; then
        log_error "PDB file not found for ${protein}: ${pdb_file}"
        return 1
    fi

    mkdir -p "$input_dir"

    if ls "${input_dir}"/*.csv >/dev/null 2>&1; then
        log_info "Using existing Vina input chunks for ${protein}: ${input_dir}"
        return 0
    fi

    log_info "Generating Vina input chunks for ${protein}"

    conda activate general
    python "$SPLIT_CSV_EXE" \
        --protein-name "$protein" \
        --pdb-file "$pdb_file" \
        --smiles-file "$selected_csv" \
        --chunk-size 100 \
        --output-dir "$input_dir"

    return 0
}

# Function to submit a single Vina docking job
submit_vina_job() {
    local protein=$1
    local part=$2
    local receptor_pdbqt=$3

    local FINE_DIR="${TASK_ROOT}/${protein}/fine_screening"
    local INPUT_CSV="${FINE_DIR}/Vina/input/${part}.csv"
    local OUTPUT_DIR="${FINE_DIR}/Vina/output/${part}"
    local TOKEN_FILE="${FINE_DIR}/Vina/output/${part}.done"

    # Check if input CSV exists
    if [ ! -f "$INPUT_CSV" ]; then
        log_error "Input CSV not found: $INPUT_CSV"
        return 1
    fi

    # Check if receptor file exists
    if [ ! -f "$receptor_pdbqt" ]; then
        log_error "Receptor PDBQT not found: $receptor_pdbqt"
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

    cat > "$JOB_SCRIPT" <<EOF_JOB
#!/bin/bash
#SBATCH --job-name=vina_${protein}_${part}
EOF_JOB

    # Add constraint if specified
    if [ -n "$CONSTRAINT" ]; then
        echo "#SBATCH --constraint=${CONSTRAINT}" >> "$JOB_SCRIPT"
    fi

    # Add GPU request if specified
    if [ -n "$GPU_REQUEST" ]; then
        echo "#SBATCH ${GPU_REQUEST}" >> "$JOB_SCRIPT"
    fi

    cat >> "$JOB_SCRIPT" <<EOF_JOB
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
conda activate ${VINA_CONDA_ENV}

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Run Vina docking
echo "Starting Vina docking..."
python "${VINA_EXE}" \
    --smiles "${INPUT_CSV}" \
    --receptor-pdbqt "${receptor_pdbqt}" \
    --boxes '${BOXES_JSON}' \
    --output "${OUTPUT_DIR}" \
    --smiles-col "ligand_description" \
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
EOF_JOB

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

log_info "Starting Vina workflow: receptor preparation -> input generation -> job submission"

if [ ! -f "$SEQS_CSV" ]; then
    log_error "sequences.csv not found: $SEQS_CSV"
    exit 1
fi

if [ ! -f "$SPLIT_CSV_EXE" ]; then
    log_error "split_csv.py not found: $SPLIT_CSV_EXE"
    exit 1
fi

if [ ! -f "$PREP_RECEPTOR_SCRIPT" ]; then
    log_error "prepare_vina_receptor.sh not found: $PREP_RECEPTOR_SCRIPT"
    exit 1
fi

# Get list of proteins to process
PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    # Step 1: Prepare receptor (required for all Vina tasks for this protein)
    if ! RECEPTOR_PDBQT=$(prepare_receptor "$PROTEIN"); then
        log_error "Receptor preparation failed for ${PROTEIN}; terminating Vina tasks for this protein"
        continue
    fi

    # Step 2: Generate Vina inputs
    if ! generate_vina_inputs "$PROTEIN"; then
        log_error "Input generation failed for ${PROTEIN}; terminating Vina tasks for this protein"
        continue
    fi

    INPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/input"

    # Find all input CSV files (parts)
    mapfile -t INPUT_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.csv" | sort)

    if [ ${#INPUT_FILES[@]} -eq 0 ]; then
        log_error "No input CSV files found in ${INPUT_DIR}; terminating Vina tasks for ${PROTEIN}"
        continue
    fi

    N_PARTS=${#INPUT_FILES[@]}
    log_info "Found ${N_PARTS} parts to process for ${PROTEIN}"

    # Clear previous job list
    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/output"
    mkdir -p "$OUTPUT_DIR"
    > "${OUTPUT_DIR}/submitted_jobs.txt"

    # Step 3: Submit docking jobs for each part
    SUBMITTED=0
    for INPUT_FILE in "${INPUT_FILES[@]}"; do
        # Extract part name (e.g., "input_part_0" from "/path/to/input_part_0.csv")
        PART=$(basename "$INPUT_FILE" .csv)

        if submit_vina_job "$PROTEIN" "$PART" "$RECEPTOR_PDBQT"; then
            SUBMITTED=$((SUBMITTED + 1))
        fi
    done

    log_info "Submitted ${SUBMITTED} docking jobs for protein ${PROTEIN}"
    log_info "Job IDs saved to: ${OUTPUT_DIR}/submitted_jobs.txt"
done

log_info "Vina workflow submission finished"
log_info "Monitor jobs with: squeue -u \$USER"
log_info "Check job outputs in: \${TASK_ROOT}/\${PROTEIN}/fine_screening/Vina/output/"
