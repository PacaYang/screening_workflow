#!/bin/bash
#
# RoseTTAFold-All-Atom Batch Screening Script (Stage 2: Protein-Ligand Prediction)
# Submits SLURM jobs for large-scale protein-ligand screening.
# Requires protein folding to be completed first (run_rosettafold_prefold.sh).
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

# RoseTTAFold-All-Atom paths
RFAA_ROOT="/home/ubuntu/Applications/RoseTTAFold-All-Atom"
RFAA_CONDA_ENV="RFAA"

# Number of batches for protein-ligand predictions
N_BATCHES=10

# SLURM configuration for protein-ligand prediction
LIGAND_TIME_LIMIT="48:00:00"
LIGAND_MEMORY="15G"
LIGAND_CPUS=4
LIGAND_GPU_REQUEST=""
LIGAND_CONSTRAINT="g5.xlarge"

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

# Function to submit protein-ligand batch job
submit_protein_ligand_batch() {
    local protein=$1
    local batch_id=$2

    local RFAA_DIR="${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold"
    local FOLD_DIR="${RFAA_DIR}/protein_folding"
    local LIGAND_DIR="${RFAA_DIR}/protein_ligand"
    local OUTPUT_DIR="${LIGAND_DIR}/output"
    local CONFIG_DIR="${LIGAND_DIR}/config"
    local LOG_DIR="${LIGAND_DIR}/logs"
    local TOKEN_DIR="${OUTPUT_DIR}/token"
    local TOKEN_FILE="${TOKEN_DIR}/batch_${batch_id}.done"
    local SELECTED_CSV="${TASK_ROOT}/${protein}/initial_screening/selected.csv"
    local FASTA_FILE="${FOLD_DIR}/input/${protein}.fasta"
    local PROTEIN_FOLD_OUTPUT="${FOLD_DIR}/output"

    # Check if already completed
    if [ -f "$TOKEN_FILE" ]; then
        log_info "Batch ${batch_id} already completed for ${protein}, skipping"
        return 0
    fi

    # Create directories
    mkdir -p "$OUTPUT_DIR" "$CONFIG_DIR" "$LOG_DIR" "$TOKEN_DIR"

    # Check if selected CSV exists
    if [ ! -f "$SELECTED_CSV" ]; then
        log_error "Selected CSV not found: $SELECTED_CSV"
        return 1
    fi

    # Create SLURM job script
    local JOB_SCRIPT="${LOG_DIR}/slurm_batch_${batch_id}.sh"

    # Create base script without dependency line
    cat > "$JOB_SCRIPT" <<'EOFSCRIPT'
#!/bin/bash
#SBATCH --job-name=rfaa_ligand_PROTEIN_PLACEHOLDER_bBATCH_ID_PLACEHOLDER
#SBATCH --time=TIME_LIMIT_PLACEHOLDER
#SBATCH --mem=MEMORY_PLACEHOLDER
#SBATCH --cpus-per-task=CPUS_PLACEHOLDER
#SBATCH GPU_REQUEST_PLACEHOLDER
#SBATCH --constraint=CONSTRAINT_PLACEHOLDER
#SBATCH --output=LOG_DIR_PLACEHOLDER/slurm_batch_BATCH_ID_PLACEHOLDER_%j.out
#SBATCH --error=LOG_DIR_PLACEHOLDER/slurm_batch_BATCH_ID_PLACEHOLDER_%j.err

set -e

echo "Job started at: $(date)"
echo "Processing protein: PROTEIN_PLACEHOLDER, batch: BATCH_ID_PLACEHOLDER"

# Activate environment
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate RFAA_CONDA_ENV_PLACEHOLDER

# Set database paths for MSA generation
export DB_UR30="/shared/programs/RFAA_data/UniRef30_2020_06/UniRef30_2020_06"
export DB_BFD="/shared/programs/RFAA_data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"

# Calculate batch range
N_COMPOUNDS=$(tail -n +2 SELECTED_CSV_PLACEHOLDER | wc -l)
BATCH_SIZE=$(( (N_COMPOUNDS + N_BATCHES_PLACEHOLDER - 1) / N_BATCHES_PLACEHOLDER ))
START_IDX=$((BATCH_ID_PLACEHOLDER * BATCH_SIZE))
END_IDX=$((START_IDX + BATCH_SIZE))
[ $END_IDX -gt $N_COMPOUNDS ] && END_IDX=$N_COMPOUNDS

echo "Processing compounds ${START_IDX} to $((END_IDX - 1)) (total: ${N_COMPOUNDS})"

# Local temporary storage
LOCAL_OUT=/tmp/rfaa_PROTEIN_PLACEHOLDER_BATCH_ID_PLACEHOLDER_${SLURM_JOB_ID}
mkdir -p $LOCAL_OUT

# Change to RFAA directory (required for relative paths)
cd RFAA_ROOT_PLACEHOLDER

# Process each compound in batch
SUCCESS_COUNT=0
FAIL_COUNT=0

for i in $(seq $START_IDX $((END_IDX - 1))); do
    echo "Processing compound index: $i"

    # Extract SMILES and compound ID for compound i
    COMPOUND_DATA=$(python3 <<PYEOF
import pandas as pd
import sys
try:
    df = pd.read_csv("SELECTED_CSV_PLACEHOLDER")
    idx = ${i}
    if idx >= len(df):
        print("ERROR: Index out of range", file=sys.stderr)
        sys.exit(1)
    row = df.iloc[idx]
    # Try common column names for SMILES
    smiles = None
    for col in ['smiles', 'SMILES', 'Smiles']:
        if col in df.columns:
            smiles = row[col]
            break
    if smiles is None:
        print("ERROR: No SMILES column found", file=sys.stderr)
        sys.exit(1)
    # Try common column names for compound ID
    compound_id = None
    for col in ['compound_id', 'id', 'name', 'ID', 'Name']:
        if col in df.columns:
            compound_id = row[col]
            break
    if compound_id is None:
        compound_id = f"compound_{idx}"
    print(f"{compound_id}|{smiles}")
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
PYEOF
)

    if [ $? -ne 0 ]; then
        echo "Failed to extract data for compound $i"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    COMPOUND_ID=$(echo "$COMPOUND_DATA" | cut -d'|' -f1)
    SMILES=$(echo "$COMPOUND_DATA" | cut -d'|' -f2)

    echo "  Compound ID: $COMPOUND_ID"
    echo "  SMILES: $SMILES"

    # Generate config for this compound
    COMPOUND_CONFIG="${LOCAL_OUT}/config_${COMPOUND_ID}.yaml"
    COMPOUND_OUTPUT="${LOCAL_OUT}/${COMPOUND_ID}"
    mkdir -p "$COMPOUND_OUTPUT"

    # Copy pre-computed MSA files from Stage 1 to avoid regenerating them
    PROTEIN_MSA_DIR="PROTEIN_FOLD_OUTPUT_PLACEHOLDER/PROTEIN_PLACEHOLDER_fold/A"
    COMPOUND_MSA_DIR="${COMPOUND_OUTPUT}/PROTEIN_PLACEHOLDER_ligand_${COMPOUND_ID}/A"
    mkdir -p "$COMPOUND_MSA_DIR"

    if [ -d "$PROTEIN_MSA_DIR" ]; then
        echo "  Copying pre-computed MSA files from Stage 1"
        cp "$PROTEIN_MSA_DIR"/t000_.msa0.a3m "$COMPOUND_MSA_DIR/" 2>/dev/null || true
        cp "$PROTEIN_MSA_DIR"/t000_.hhr "$COMPOUND_MSA_DIR/" 2>/dev/null || true
        cp "$PROTEIN_MSA_DIR"/t000_.atab "$COMPOUND_MSA_DIR/" 2>/dev/null || true
        cp "$PROTEIN_MSA_DIR"/t000_.ss2 "$COMPOUND_MSA_DIR/" 2>/dev/null || true
    else
        echo "  Warning: Pre-computed MSA directory not found: $PROTEIN_MSA_DIR"
    fi

    cat > "$COMPOUND_CONFIG" <<EOF
defaults:
  - base

job_name: "PROTEIN_PLACEHOLDER_ligand_${COMPOUND_ID}"
output_path: "${COMPOUND_OUTPUT}"

protein_inputs:
  A:
    fasta_file: "FASTA_FILE_PLACEHOLDER"

sm_inputs:
  B:
    input: "${SMILES}"
    input_type: "smiles"
EOF

    # Run inference for this compound
    if python -m rf2aa.run_inference \
        --config-dir "${LOCAL_OUT}" \
        --config-name "config_${COMPOUND_ID}" 2>&1 | tee "${LOCAL_OUT}/${COMPOUND_ID}.log"; then
        echo "  SUCCESS: Compound $COMPOUND_ID"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "  FAILED: Compound $COMPOUND_ID"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
done

echo "Batch processing complete:"
echo "  Success: ${SUCCESS_COUNT}"
echo "  Failed: ${FAIL_COUNT}"

# Compress and copy results
if [ $SUCCESS_COUNT -gt 0 ]; then
    tar -czf OUTPUT_DIR_PLACEHOLDER/batch_BATCH_ID_PLACEHOLDER.tar.gz -C $LOCAL_OUT .
    echo "Results compressed to OUTPUT_DIR_PLACEHOLDER/batch_BATCH_ID_PLACEHOLDER.tar.gz"
fi

# Create completion token
touch TOKEN_FILE_PLACEHOLDER

# Cleanup
rm -rf $LOCAL_OUT

echo "Job completed at: $(date)"
EOFSCRIPT

    # Replace placeholders
    sed -i "s|PROTEIN_PLACEHOLDER|${protein}|g" "$JOB_SCRIPT"
    sed -i "s|BATCH_ID_PLACEHOLDER|${batch_id}|g" "$JOB_SCRIPT"
    sed -i "s|TIME_LIMIT_PLACEHOLDER|${LIGAND_TIME_LIMIT}|g" "$JOB_SCRIPT"
    sed -i "s|MEMORY_PLACEHOLDER|${LIGAND_MEMORY}|g" "$JOB_SCRIPT"
    sed -i "s|CPUS_PLACEHOLDER|${LIGAND_CPUS}|g" "$JOB_SCRIPT"
    sed -i "s|GPU_REQUEST_PLACEHOLDER|${LIGAND_GPU_REQUEST}|g" "$JOB_SCRIPT"
    sed -i "s|CONSTRAINT_PLACEHOLDER|${LIGAND_CONSTRAINT}|g" "$JOB_SCRIPT"
    sed -i "s|LOG_DIR_PLACEHOLDER|${LOG_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|RFAA_CONDA_ENV_PLACEHOLDER|${RFAA_CONDA_ENV}|g" "$JOB_SCRIPT"
    sed -i "s|RFAA_ROOT_PLACEHOLDER|${RFAA_ROOT}|g" "$JOB_SCRIPT"
    sed -i "s|SELECTED_CSV_PLACEHOLDER|${SELECTED_CSV}|g" "$JOB_SCRIPT"
    sed -i "s|N_BATCHES_PLACEHOLDER|${N_BATCHES}|g" "$JOB_SCRIPT"
    sed -i "s|FASTA_FILE_PLACEHOLDER|${FASTA_FILE}|g" "$JOB_SCRIPT"
    sed -i "s|OUTPUT_DIR_PLACEHOLDER|${OUTPUT_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|TOKEN_FILE_PLACEHOLDER|${TOKEN_FILE}|g" "$JOB_SCRIPT"
    sed -i "s|PROTEIN_FOLD_OUTPUT_PLACEHOLDER|${PROTEIN_FOLD_OUTPUT}|g" "$JOB_SCRIPT"

    # Submit job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted protein-ligand batch ${batch_id} for ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${RFAA_DIR}/submitted_jobs.txt"
        return 0
    else
        log_error "Failed to submit protein-ligand batch ${batch_id} for ${protein}"
        return 1
    fi
}

# ============================================================================
# Main Script
# ============================================================================

log_info "Starting RoseTTAFold-All-Atom protein-ligand batch screening"
log_info "Task root: ${TASK_ROOT}"
log_info "Number of batches: ${N_BATCHES}"

# Pre-flight checks
if [ ! -d "$RFAA_ROOT" ]; then
    log_error "RoseTTAFold-All-Atom directory not found: ${RFAA_ROOT}"
    exit 1
fi

# Get protein list
PROTEINS=$(get_proteins)
log_info "Processing proteins: ${PROTEINS}"

# Process each protein
for PROTEIN in $PROTEINS; do
    log_info "=========================================="
    log_info "Processing protein: ${PROTEIN}"
    log_info "=========================================="

    # Check if protein folding is completed (prerequisite)
    FOLD_TOKEN="${TASK_ROOT}/${PROTEIN}/fine_screening/RoseTTAFold/protein_folding/output/protein_fold.done"
    if [ ! -f "$FOLD_TOKEN" ]; then
        log_error "Protein folding not completed for ${PROTEIN}: ${FOLD_TOKEN}"
        log_error "Please run run_rosettafold_prefold.sh first"
        continue
    fi

    # Validate protein directory structure
    PROTEIN_DIR="${TASK_ROOT}/${PROTEIN}"
    SELECTED_CSV="${PROTEIN_DIR}/initial_screening/selected.csv"

    if [ ! -f "$SELECTED_CSV" ]; then
        log_error "Selected compounds CSV not found: ${SELECTED_CSV}"
        continue
    fi

    # Count compounds
    N_COMPOUNDS=$(tail -n +2 "$SELECTED_CSV" | wc -l)
    log_info "Found ${N_COMPOUNDS} compounds for ${PROTEIN}"

    if [ "$N_COMPOUNDS" -eq 0 ]; then
        log_error "No compounds found in ${SELECTED_CSV}"
        continue
    fi

    BATCH_SIZE=$(( (N_COMPOUNDS + N_BATCHES - 1) / N_BATCHES ))
    log_info "Batch size: ${BATCH_SIZE} compounds per batch"

    RFAA_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/RoseTTAFold"
    > "${RFAA_DIR}/submitted_jobs.txt"

    SUBMITTED=0
    for batch_id in $(seq 0 $((N_BATCHES - 1))); do
        START_IDX=$((batch_id * BATCH_SIZE))
        if [ $START_IDX -ge $N_COMPOUNDS ]; then
            break
        fi

        if submit_protein_ligand_batch "$PROTEIN" "$batch_id"; then
            SUBMITTED=$((SUBMITTED + 1))
        fi
    done

    log_info "Submitted ${SUBMITTED} batch jobs for ${PROTEIN}"
done

log_info "=========================================="
log_info "All jobs submitted"
log_info "=========================================="
log_info "Monitor job status with: squeue -u \$USER"
log_info "Check logs in: \${TASK_ROOT}/\${PROTEIN}/fine_screening/RoseTTAFold/protein_ligand/logs/"

exit 0
