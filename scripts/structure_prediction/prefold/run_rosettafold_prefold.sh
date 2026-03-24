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
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

# Load configuration from environment or use defaults
TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"

# RoseTTAFold-All-Atom paths
RFAA_ROOT="/home/yangl_pacagen_com/Applications/RoseTTAFold-All-Atom"
RFAA_BASE_CONFIG="${RFAA_ROOT}/rf2aa/config/inference/base.yaml"
RFAA_WEIGHTS="${MASTER_RFAA_WEIGHTS:-/home/yangl_pacagen_com/Applications/model_weights/RoseTTAFold/RFAA_paper_weights.pt}"
ROSETTA_DB_UR30="${MASTER_ROSETTA_DB_UR30:-/home/yangl_pacagen_com/Applications/model_weights/rosetta_db/UniRef30_2020_06/UniRef30_2020_06}"
ROSETTA_DB_BFD="${MASTER_ROSETTA_DB_BFD:-/home/yangl_pacagen_com/Applications/model_weights/rosetta_db/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt}"
RFAA_CONDA_ENV="${MASTER_RFAA_CONDA_ENV:-RFAA}"
PREFOLD_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAFE_MSA_WRAPPER="${PREFOLD_SCRIPT_DIR}/make_msa_safe.sh"
SAFE_MSA_COMMAND_REL=""

# SLURM configuration for protein folding
PROTEIN_TIME_LIMIT="${PROTEIN_TIME_LIMIT:-48:00:00}"
PROTEIN_MEMORY="${PROTEIN_MEMORY:-64G}"
PROTEIN_CPUS=${PROTEIN_CPUS:-16}
PROTEIN_EXCLUSIVE="${PROTEIN_EXCLUSIVE:-1}"
PROTEIN_GPU_REQUEST="${PROTEIN_GPU_REQUEST:---gres=gpu:1}"
PROTEIN_PARTITION="${PROTEIN_PARTITION:-g232}"

# ============================================================================
# Functions
# ============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*" >&2
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN: $*" >&2
}

resolve_safe_msa_command_relpath() {
    local wrapper_path=$1
    realpath --relative-to "$RFAA_ROOT" "$wrapper_path"
}

# Convert SLURM memory notation (e.g., 64G, 64000M, 15360) to integer GB.
parse_memory_to_gb() {
    local memory=$1
    python3 - "$memory" <<'PYEOF'
import math
import re
import sys

raw = sys.argv[1].strip().upper()
match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s*([KMGT]?)B?", raw)
if not match:
    sys.exit(1)

value = float(match.group(1))
unit = match.group(2)

# Bare numeric values in Slurm are megabytes.
factor_to_gb = {
    "": 1.0 / 1024.0,
    "K": 1.0 / (1024.0 * 1024.0),
    "M": 1.0 / 1024.0,
    "G": 1.0,
    "T": 1024.0,
}

print(int(math.ceil(value * factor_to_gb[unit])))
PYEOF
}

# Read the RoseTTAFold base config's MSA memory target (GB).
detect_rfaa_msa_mem_gb() {
    [ -f "$RFAA_BASE_CONFIG" ] || return 1
    python3 - "$RFAA_BASE_CONFIG" <<'PYEOF'
import re
import sys

config_path = sys.argv[1]
with open(config_path, "r", encoding="utf-8") as handle:
    for line in handle:
        if re.match(r"^\s*mem\s*:", line):
            value = line.split(":", 1)[1].split("#", 1)[0].strip().strip('"').strip("'")
            if not value:
                break
            try:
                print(int(float(value)))
                sys.exit(0)
            except ValueError:
                break
sys.exit(1)
PYEOF
}

# Ensure SLURM allocation cannot be lower than the MSA memory target.
align_prefold_memory_with_rfaa_config() {
    local rfaa_msa_mem_gb
    local slurm_mem_gb

    if ! rfaa_msa_mem_gb="$(detect_rfaa_msa_mem_gb)"; then
        log_warn "Could not detect MSA mem from ${RFAA_BASE_CONFIG}; using PROTEIN_MEMORY=${PROTEIN_MEMORY}"
        return 0
    fi

    if ! slurm_mem_gb="$(parse_memory_to_gb "${PROTEIN_MEMORY}")"; then
        log_warn "Could not parse PROTEIN_MEMORY='${PROTEIN_MEMORY}'; expected format like 64G or 64000M"
        return 0
    fi

    if [ "${slurm_mem_gb}" -lt "${rfaa_msa_mem_gb}" ]; then
        log_warn "PROTEIN_MEMORY=${PROTEIN_MEMORY} is below RoseTTAFold MSA mem (${rfaa_msa_mem_gb}G); bumping to ${rfaa_msa_mem_gb}G"
        PROTEIN_MEMORY="${rfaa_msa_mem_gb}G"
    fi

    log_info "Using PROTEIN_MEMORY=${PROTEIN_MEMORY} (RoseTTAFold MSA mem=${rfaa_msa_mem_gb}G)"
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
    protein_row = df[df['protein_name'] == "${protein}"]

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
    local msa_command=$5

    log_info "Generating protein folding config for ${protein}"

    local pdb100_db="/home/yangl_pacagen_com/Applications/model_weights/rosetta_db/pdb100_2021Mar03/pdb100_2021Mar03"

    cat > "$config_file" <<EOF
defaults:
  - base

job_name: "${protein}_fold"
output_path: "${output_dir}"
checkpoint_path: "${RFAA_WEIGHTS}"

database_params:
  command: "${msa_command}"
  hhdb: "${pdb100_db}"

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
    generate_protein_fold_config "$protein" "$CONFIG_FILE" "$FASTA_FILE" "$OUTPUT_DIR" "$SAFE_MSA_COMMAND_REL"

    # Create SLURM job script
    local JOB_SCRIPT="${LOG_DIR}/slurm_protein_fold.sh"

    cat > "$JOB_SCRIPT" <<'EOFSCRIPT'
#!/bin/bash
#SBATCH --job-name=rfaa_fold_PROTEIN_PLACEHOLDER
#SBATCH --time=TIME_LIMIT_PLACEHOLDER
#SBATCH --mem=MEMORY_PLACEHOLDER
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=CPUS_PLACEHOLDER
#SBATCH EXCLUSIVE_PLACEHOLDER
#SBATCH --partition=PARTITION_PLACEHOLDER
#SBATCH --output=LOG_DIR_PLACEHOLDER/slurm_protein_%j.out
#SBATCH --error=LOG_DIR_PLACEHOLDER/slurm_protein_%j.err

set -e

echo "Job started at: $(date)"
echo "Folding protein: PROTEIN_PLACEHOLDER"

# Set database paths
export DB_UR30="DB_UR30_PLACEHOLDER"
export DB_BFD="DB_BFD_PLACEHOLDER"
export RFAA_ROOT="RFAA_ROOT_PLACEHOLDER"

# Activate environment
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate RFAA_CONDA_ENV_PLACEHOLDER

# Change to RFAA directory (required for relative paths)
cd RFAA_ROOT_PLACEHOLDER

# Run protein folding
echo "Config path: CONFIG_DIR_PLACEHOLDER/protein_fold.yaml"
echo "MSA command override: MSA_COMMAND_PLACEHOLDER"
if ! python -m rf2aa.run_inference \
    --config-dir CONFIG_DIR_PLACEHOLDER \
    --config-name protein_fold; then
    echo "ERROR: rf2aa.run_inference failed for PROTEIN_PLACEHOLDER"
    exit 1
fi

# Validate required prefold artifacts before creating completion tokens
REQUIRED_OUTPUTS=(
    "OUTPUT_DIR_PLACEHOLDER/PROTEIN_PLACEHOLDER_fold.pdb"
    "OUTPUT_DIR_PLACEHOLDER/PROTEIN_PLACEHOLDER_fold_aux.pt"
    "OUTPUT_DIR_PLACEHOLDER/PROTEIN_PLACEHOLDER_fold/A/t000_.msa0.a3m"
    "OUTPUT_DIR_PLACEHOLDER/PROTEIN_PLACEHOLDER_fold/A/t000_.ss2"
    "OUTPUT_DIR_PLACEHOLDER/PROTEIN_PLACEHOLDER_fold/A/t000_.hhr"
    "OUTPUT_DIR_PLACEHOLDER/PROTEIN_PLACEHOLDER_fold/A/t000_.atab"
)
for required_path in "${REQUIRED_OUTPUTS[@]}"; do
    if [ ! -s "$required_path" ]; then
        echo "ERROR: Missing expected prefold output: $required_path"
        exit 1
    fi
done

# Create completion tokens only after full validation
touch TOKEN_FILE_PLACEHOLDER

echo "Job completed at: $(date)"
EOFSCRIPT

    # Replace placeholders
    sed -i "s|PROTEIN_PLACEHOLDER|${protein}|g" "$JOB_SCRIPT"
    sed -i "s|TIME_LIMIT_PLACEHOLDER|${PROTEIN_TIME_LIMIT}|g" "$JOB_SCRIPT"
    sed -i "s|MEMORY_PLACEHOLDER|${PROTEIN_MEMORY}|g" "$JOB_SCRIPT"
    sed -i "s|CPUS_PLACEHOLDER|${PROTEIN_CPUS}|g" "$JOB_SCRIPT"
    sed -i "s|GPU_REQUEST_PLACEHOLDER|${PROTEIN_GPU_REQUEST}|g" "$JOB_SCRIPT"
    # Set --exclusive if explicitly requested, otherwise remove the placeholder line
    if [ "${PROTEIN_EXCLUSIVE}" = "1" ]; then
        sed -i "s|#SBATCH EXCLUSIVE_PLACEHOLDER|#SBATCH --exclusive|g" "$JOB_SCRIPT"
    else
        sed -i "/#SBATCH EXCLUSIVE_PLACEHOLDER/d" "$JOB_SCRIPT"
    fi
    sed -i "s|LOG_DIR_PLACEHOLDER|${LOG_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|RFAA_CONDA_ENV_PLACEHOLDER|${RFAA_CONDA_ENV}|g" "$JOB_SCRIPT"
    sed -i "s|RFAA_ROOT_PLACEHOLDER|${RFAA_ROOT}|g" "$JOB_SCRIPT"
    sed -i "s|CONFIG_DIR_PLACEHOLDER|${CONFIG_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|OUTPUT_DIR_PLACEHOLDER|${OUTPUT_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|TOKEN_FILE_PLACEHOLDER|${TOKEN_FILE}|g" "$JOB_SCRIPT"
    sed -i "s|DB_UR30_PLACEHOLDER|${ROSETTA_DB_UR30}|g" "$JOB_SCRIPT"
    sed -i "s|DB_BFD_PLACEHOLDER|${ROSETTA_DB_BFD}|g" "$JOB_SCRIPT"
    sed -i "s|PARTITION_PLACEHOLDER|${PROTEIN_PARTITION}|g" "$JOB_SCRIPT"
    sed -i "s|MSA_COMMAND_PLACEHOLDER|${SAFE_MSA_COMMAND_REL}|g" "$JOB_SCRIPT"

    # Submit job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted protein folding job for ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${RFAA_DIR}/job_ids.txt"
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

if [ ! -x "$SAFE_MSA_WRAPPER" ]; then
    log_error "Safe MSA wrapper not found or not executable: ${SAFE_MSA_WRAPPER}"
    exit 1
fi

if ! SAFE_MSA_COMMAND_REL="$(resolve_safe_msa_command_relpath "$SAFE_MSA_WRAPPER")"; then
    log_error "Failed to resolve safe MSA wrapper relative to RFAA root: ${SAFE_MSA_WRAPPER}"
    exit 1
fi

if [[ "$SAFE_MSA_COMMAND_REL" = /* ]]; then
    log_error "Safe MSA command must be a relative path from ${RFAA_ROOT}, got: ${SAFE_MSA_COMMAND_REL}"
    exit 1
fi

log_info "Using safe MSA command override: ${SAFE_MSA_COMMAND_REL}"

# Ensure memory requested from SLURM is compatible with RoseTTAFold MSA settings.
align_prefold_memory_with_rfaa_config

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
