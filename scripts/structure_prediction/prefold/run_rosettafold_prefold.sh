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

# Hard memory cap for the MSA step. hhblits/hhsearch are invoked with -maxmem
# equal to RFAA database_params.mem; that value MUST stay strictly below the SLURM
# allocation or the cgroup OOM-kills the job (hhblits routinely overshoots -maxmem).
# We pass (SLURM_mem_GB - MSA_HEADROOM_GB) as RFAA's mem, floored at MSA_MIN_MEM_GB.
MSA_HEADROOM_GB="${MSA_HEADROOM_GB:-16}"
MSA_MIN_MEM_GB="${MSA_MIN_MEM_GB:-8}"
# Resolved by enforce_msa_memory_headroom(); fed into each fold config.
RFAA_MSA_MEM_GB=""
# Path to shared chain-parsing helpers (seq_utils.parse_chains / chain_ids).
SEQ_UTILS_DIR="$(cd "${PREFOLD_SCRIPT_DIR}/../../input" && pwd)"

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

# Enforce a hard memory headroom between the MSA tool's -maxmem target and the SLURM
# allocation. hhblits/hhsearch take -maxmem from RFAA database_params.mem; if that
# equals the SLURM --mem, the cgroup OOM-kills the overshoot. We instead set RFAA's
# mem to (SLURM_mem_GB - MSA_HEADROOM_GB), floored at MSA_MIN_MEM_GB, and feed that
# into every fold config's database_params.mem. SLURM --mem stays the outer hard wall.
enforce_msa_memory_headroom() {
    local slurm_mem_gb
    local capped

    if ! slurm_mem_gb="$(parse_memory_to_gb "${PROTEIN_MEMORY}")"; then
        log_warn "Could not parse PROTEIN_MEMORY='${PROTEIN_MEMORY}'; expected format like 64G or 64000M"
        log_warn "Falling back to MSA mem=${MSA_MIN_MEM_GB}G"
        RFAA_MSA_MEM_GB="${MSA_MIN_MEM_GB}"
        return 0
    fi

    capped=$(( slurm_mem_gb - MSA_HEADROOM_GB ))
    if [ "${capped}" -lt "${MSA_MIN_MEM_GB}" ]; then
        capped="${MSA_MIN_MEM_GB}"
        # If even the floor leaves no headroom, bump the SLURM allocation so the
        # invariant (maxmem < slurm_mem) still holds.
        if [ "$(( capped + MSA_HEADROOM_GB ))" -gt "${slurm_mem_gb}" ]; then
            local bumped=$(( capped + MSA_HEADROOM_GB ))
            log_warn "PROTEIN_MEMORY=${PROTEIN_MEMORY} too small for MSA floor ${MSA_MIN_MEM_GB}G + headroom ${MSA_HEADROOM_GB}G; bumping SLURM mem to ${bumped}G"
            PROTEIN_MEMORY="${bumped}G"
            slurm_mem_gb="${bumped}"
        fi
    fi

    RFAA_MSA_MEM_GB="${capped}"
    log_info "MSA memory cap: slurm_mem=${slurm_mem_gb}G maxmem=${RFAA_MSA_MEM_GB}G headroom=${MSA_HEADROOM_GB}G"
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
# Write one FASTA per protein chain into input_dir, named <protein>_<CID>.fasta
# (A, B, C...). The comma-separated `sequence` cell is split via seq_utils.parse_chains
# (B2/MCP convention). Echoes the chain count on stdout so the caller can iterate.
extract_protein_fasta() {
    local protein=$1
    local input_dir=$2

    log_info "Extracting per-chain FASTA sequences for ${protein}"

    PYTHONPATH="${SEQ_UTILS_DIR}:${PYTHONPATH:-}" python3 - "$SEQS_CSV" "$protein" "$input_dir" <<'PYEOF'
import sys
import pandas as pd
from seq_utils import parse_chains, chain_ids

seqs_csv, protein, input_dir = sys.argv[1], sys.argv[2], sys.argv[3]
try:
    df = pd.read_csv(seqs_csv)
    # Target name column: prefer 'name' (B2/MCP), fall back to legacy 'protein_name'.
    name_col = next((c for c in ['name', 'protein_name'] if c in df.columns), None)
    if name_col is None:
        print("ERROR: No name/protein_name column in sequences.csv", file=sys.stderr)
        sys.exit(1)
    row = df[df[name_col] == protein]
    if row.empty:
        print(f"ERROR: Protein {protein} not found in sequences.csv", file=sys.stderr)
        sys.exit(1)
    seq_col = next((c for c in ['sequence', 'fasta', 'seq', 'protein_sequence'] if c in df.columns), None)
    if seq_col is None:
        print("ERROR: No sequence column found in sequences.csv", file=sys.stderr)
        sys.exit(1)
    chains = parse_chains(row[seq_col].iloc[0])
    if not chains:
        print(f"ERROR: No chain sequences parsed for {protein}", file=sys.stderr)
        sys.exit(1)
    for cid, seq in zip(chain_ids(len(chains)), chains):
        with open(f"{input_dir}/{protein}_{cid}.fasta", 'w') as f:
            f.write(f">{protein}_{cid}\n{seq}\n")
    # Print chain count on the last line for the bash caller.
    print(len(chains))
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
PYEOF
}

# Function to generate Hydra config for protein folding (multi-chain).
generate_protein_fold_config() {
    local protein=$1
    local config_file=$2
    local input_dir=$3
    local output_dir=$4
    local msa_command=$5
    local n_chains=$6

    log_info "Generating protein folding config for ${protein} (${n_chains} chains)"

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
  mem: ${RFAA_MSA_MEM_GB}

protein_inputs:
EOF

    # One protein_inputs entry per chain (A, B, C...), each pointing at its FASTA.
    # make_msa writes each chain's MSA to output/${protein}_fold/<CID>/t000_.*
    PYTHONPATH="${SEQ_UTILS_DIR}:${PYTHONPATH:-}" python3 - "$protein" "$input_dir" "$n_chains" >> "$config_file" <<'PYEOF'
import sys
from seq_utils import chain_ids
protein, input_dir, n_chains = sys.argv[1], sys.argv[2], int(sys.argv[3])
for cid in chain_ids(n_chains):
    print(f"  {cid}:")
    print(f"    fasta_file: \"{input_dir}/{protein}_{cid}.fasta\"")
PYEOF

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
    local CONFIG_FILE="${CONFIG_DIR}/protein_fold.yaml"

    # Check if already completed
    if [ -f "$TOKEN_FILE" ]; then
        log_info "Protein folding already completed for ${protein}, skipping"
        return 0
    fi

    # Create directories
    mkdir -p "$INPUT_DIR" "$OUTPUT_DIR" "$CONFIG_DIR" "$LOG_DIR"

    # Write one FASTA per chain; capture the chain count (last line of stdout).
    local FASTA_OUT
    if ! FASTA_OUT="$(extract_protein_fasta "$protein" "$INPUT_DIR")"; then
        log_error "Failed to extract FASTA for ${protein}: ${FASTA_OUT}"
        return 1
    fi
    local N_CHAINS
    N_CHAINS="$(echo "$FASTA_OUT" | tail -n1)"
    if ! [[ "$N_CHAINS" =~ ^[0-9]+$ ]] || [ "$N_CHAINS" -lt 1 ]; then
        log_error "Invalid chain count for ${protein}: '${N_CHAINS}'"
        return 1
    fi
    log_info "Folding ${protein} with ${N_CHAINS} protein chain(s)"

    # Generate config file (multi-chain protein_inputs)
    generate_protein_fold_config "$protein" "$CONFIG_FILE" "$INPUT_DIR" "$OUTPUT_DIR" "$SAFE_MSA_COMMAND_REL" "$N_CHAINS"

    # Build the per-chain MSA validation block injected into the job script.
    local REQUIRED_CHAIN_OUTPUTS
    REQUIRED_CHAIN_OUTPUTS="$(PYTHONPATH="${SEQ_UTILS_DIR}:${PYTHONPATH:-}" python3 - "$N_CHAINS" <<'PYEOF'
import sys
from seq_utils import chain_ids
n = int(sys.argv[1])
for cid in chain_ids(n):
    for ext in ("msa0.a3m", "ss2", "hhr", "atab"):
        print(f'    "OUTPUT_DIR_PLACEHOLDER/PROTEIN_PLACEHOLDER_fold/{cid}/t000_.{ext}"')
PYEOF
)"

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
REQUIRED_CHAIN_OUTPUTS_PLACEHOLDER
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
    # Inject the per-chain validation lines first (multi-line), using a Python helper
    # so newlines and slashes in paths don't break sed. The injected lines still
    # contain PROTEIN_/OUTPUT_DIR placeholders, substituted by the sed calls below.
    REQUIRED_CHAIN_OUTPUTS="$REQUIRED_CHAIN_OUTPUTS" python3 - "$JOB_SCRIPT" <<'PYEOF'
import os, sys
path = sys.argv[1]
block = os.environ.get("REQUIRED_CHAIN_OUTPUTS", "")
with open(path) as f:
    text = f.read()
text = text.replace("REQUIRED_CHAIN_OUTPUTS_PLACEHOLDER", block)
with open(path, "w") as f:
    f.write(text)
PYEOF
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

# Cap the MSA tool's memory below the SLURM allocation to prevent OOM kills.
enforce_msa_memory_headroom

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
