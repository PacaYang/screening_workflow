#!/bin/bash
#
# RoseTTAFold-All-Atom Protein Folding Script (Stage 1)
# Submits SLURM jobs to fold protein structures from FASTA sequences.
# Run this before run_rosettafold_batch.sh (Stage 2: protein-ligand prediction).
#

set -e

# ============================================================================
# Configuration - Update these paths as needed
# ============================================================================

# Source conda configuration
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

# Load configuration from environment or use defaults
TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"

# RoseTTAFold-All-Atom paths
RFAA_ROOT="/home/ubuntu/Applications/RoseTTAFold-All-Atom"
RFAA_WEIGHTS="${RFAA_ROOT}/RFAA_paper_weights.pt"
RFAA_CONDA_ENV="RFAA"

# SLURM configuration for protein folding
PROTEIN_TIME_LIMIT="48:00:00"
PROTEIN_MEMORY="32G"
PROTEIN_CPUS=8
PROTEIN_GPU_REQUEST="--gres=gpu:a10g:1"

# ============================================================================
# Functions
# ============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*" >&2
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

# Function to extract protein FASTA sequence from sequences.csv
extract_protein_fasta() {
    local protein=$1
    local output_fasta=$2

    log_info "Extracting FASTA sequence for ${protein}"

    python3 <<EOF
import pandas as pd
import sys

try:
    df = pd.read_csv("${SEQS_CSV}")

    # Find the protein row
    protein_row = df[df['name'] == "${protein}"]

    if protein_row.empty:
        print("ERROR: Protein ${protein} not found in sequences.csv", file=sys.stderr)
        sys.exit(1)

    # Extract sequence (try common column names)
    sequence = None
    for col in ['sequence', 'fasta', 'seq', 'protein_sequence']:
        if col in df.columns:
            sequence = protein_row[col].iloc[0]
            break

    if sequence is None:
        print("ERROR: No sequence column found in sequences.csv", file=sys.stderr)
        sys.exit(1)

    # Write FASTA file
    with open("${output_fasta}", 'w') as f:
        f.write(f">${protein}\\n")
        f.write(f"{sequence}\\n")

    print(f"FASTA written to ${output_fasta}")

except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    return $?
}

# Function to generate Hydra config for protein folding
generate_protein_fold_config() {
    local protein=$1
    local config_file=$2
    local fasta_file=$3
    local output_dir=$4

    log_info "Generating protein folding config for ${protein}"

    cat > "$config_file" <<EOF
defaults:
  - base

job_name: "${protein}_fold"
output_path: "${output_dir}"

protein_inputs:
  A:
    fasta_file: "${fasta_file}"
EOF

    log_info "Config written to ${config_file}"
}

# Function to submit protein folding job
submit_protein_fold_job() {
    local protein=$1

    local RFAA_DIR="${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold"
    local FOLD_DIR="${RFAA_DIR}/protein_folding"
    local INPUT_DIR="${FOLD_DIR}/input"
    local OUTPUT_DIR="${FOLD_DIR}/output"
    local CONFIG_DIR="${FOLD_DIR}/config"
    local LOG_DIR="${FOLD_DIR}/logs"
    local TOKEN_FILE="${OUTPUT_DIR}/protein_fold.done"
    local FASTA_FILE="${INPUT_DIR}/${protein}.fasta"
    local CONFIG_FILE="${CONFIG_DIR}/protein_fold.yaml"

    # Check if already completed
    if [ -f "$TOKEN_FILE" ]; then
        log_info "Protein folding already completed for ${protein}, skipping"
        return 0
    fi

    # Create directories
    mkdir -p "$INPUT_DIR" "$OUTPUT_DIR" "$CONFIG_DIR" "$LOG_DIR"

    # Extract FASTA sequence
    if ! extract_protein_fasta "$protein" "$FASTA_FILE"; then
        log_error "Failed to extract FASTA for ${protein}"
        return 1
    fi

    # Generate config file
    generate_protein_fold_config "$protein" "$CONFIG_FILE" "$FASTA_FILE" "$OUTPUT_DIR"

    # Create SLURM job script
    local JOB_SCRIPT="${LOG_DIR}/slurm_protein_fold.sh"

    cat > "$JOB_SCRIPT" <<'EOFSCRIPT'
#!/bin/bash
#SBATCH --job-name=rfaa_fold_PROTEIN_PLACEHOLDER
#SBATCH --time=TIME_LIMIT_PLACEHOLDER
#SBATCH --mem=MEMORY_PLACEHOLDER
#SBATCH --cpus-per-task=CPUS_PLACEHOLDER
#SBATCH GPU_REQUEST_PLACEHOLDER
#SBATCH --output=LOG_DIR_PLACEHOLDER/slurm_protein_%j.out
#SBATCH --error=LOG_DIR_PLACEHOLDER/slurm_protein_%j.err

set -e

echo "Job started at: $(date)"
echo "Folding protein: PROTEIN_PLACEHOLDER"

# Set database paths
export DB_UR30="/shared/programs/RFAA_data/UniRef30_2020_06/UniRef30_2020_06"
export DB_BFD="/shared/programs/RFAA_data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"

# Activate environment
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate RFAA_CONDA_ENV_PLACEHOLDER

# Change to RFAA directory (required for relative paths)
cd RFAA_ROOT_PLACEHOLDER

# Run protein folding
python -m rf2aa.run_inference \
    --config-dir CONFIG_DIR_PLACEHOLDER \
    --config-name protein_fold

# Create completion token
touch TOKEN_FILE_PLACEHOLDER

echo "Job completed at: $(date)"
EOFSCRIPT

    # Replace placeholders
    sed -i "s|PROTEIN_PLACEHOLDER|${protein}|g" "$JOB_SCRIPT"
    sed -i "s|TIME_LIMIT_PLACEHOLDER|${PROTEIN_TIME_LIMIT}|g" "$JOB_SCRIPT"
    sed -i "s|MEMORY_PLACEHOLDER|${PROTEIN_MEMORY}|g" "$JOB_SCRIPT"
    sed -i "s|CPUS_PLACEHOLDER|${PROTEIN_CPUS}|g" "$JOB_SCRIPT"
    sed -i "s|GPU_REQUEST_PLACEHOLDER|${PROTEIN_GPU_REQUEST}|g" "$JOB_SCRIPT"
    sed -i "s|LOG_DIR_PLACEHOLDER|${LOG_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|RFAA_CONDA_ENV_PLACEHOLDER|${RFAA_CONDA_ENV}|g" "$JOB_SCRIPT"
    sed -i "s|RFAA_ROOT_PLACEHOLDER|${RFAA_ROOT}|g" "$JOB_SCRIPT"
    sed -i "s|CONFIG_DIR_PLACEHOLDER|${CONFIG_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|TOKEN_FILE_PLACEHOLDER|${TOKEN_FILE}|g" "$JOB_SCRIPT"

    # Submit job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted protein folding job for ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${RFAA_DIR}/submitted_jobs.txt"
        return 0
    else
        log_error "Failed to submit protein folding job for ${protein}"
        return 1
    fi
}

# ============================================================================
# Main Script
# ============================================================================

log_info "Starting RoseTTAFold-All-Atom protein folding"
log_info "Task root: ${TASK_ROOT}"
log_info "RoseTTAFold root: ${RFAA_ROOT}"

# Pre-flight checks
if [ ! -d "$RFAA_ROOT" ]; then
    log_error "RoseTTAFold-All-Atom directory not found: ${RFAA_ROOT}"
    exit 1
fi

if [ ! -f "$RFAA_WEIGHTS" ]; then
    log_error "RoseTTAFold-All-Atom weights not found: ${RFAA_WEIGHTS}"
    exit 1
fi

if [ ! -f "$SEQS_CSV" ]; then
    log_error "Sequences CSV not found: ${SEQS_CSV}"
    exit 1
fi

# Check conda environment
if ! conda env list | grep -q "^${RFAA_CONDA_ENV} "; then
    log_error "Conda environment not found: ${RFAA_CONDA_ENV}"
    exit 1
fi

# Get protein list
PROTEINS=$(get_proteins)
log_info "Processing proteins: ${PROTEINS}"

SUBMITTED=0

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: ${PROTEIN}"

    if submit_protein_fold_job "$PROTEIN"; then
        SUBMITTED=$((SUBMITTED + 1))
    fi
done

log_info "Submitted ${SUBMITTED} protein folding jobs"
log_info "Monitor jobs with: squeue -u \$USER"
