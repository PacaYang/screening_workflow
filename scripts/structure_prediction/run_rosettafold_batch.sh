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
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

# Load configuration from environment or use defaults
TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
SCRIPT_ROOT="/home/yangl_pacagen_com/screening_workflow/scripts"

# RoseTTAFold-All-Atom paths
RFAA_ROOT="/home/yangl_pacagen_com/Applications/RoseTTAFold-All-Atom"
RFAA_CONDA_ENV="${MASTER_RFAA_CONDA_ENV:-RFAA}"
ROSETTA_DB_UR30="${MASTER_ROSETTA_DB_UR30:-/home/yangl_pacagen_com/Applications/model_weights/rosetta_db/UniRef30_2020_06/UniRef30_2020_06}"
ROSETTA_DB_BFD="${MASTER_ROSETTA_DB_BFD:-/home/yangl_pacagen_com/Applications/model_weights/rosetta_db/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt}"

# Number of batches for protein-ligand predictions
N_BATCHES=10

# SLURM configuration for protein-ligand prediction
LIGAND_TIME_LIMIT="48:00:00"
LIGAND_MEMORY="15G"
LIGAND_CPUS=2
LIGAND_GPU_REQUEST=""
LIGAND_PARTITION="g24"

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
    local LEDGER_FILE="${OUTPUT_DIR}/job_results.tsv"
    local LEDGER_HELPER="${SCRIPT_ROOT}/pipeline/lib/job_result_ledger.sh"
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
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=CPUS_PLACEHOLDER
#SBATCH GPU_REQUEST_PLACEHOLDER
#SBATCH --partition=PARTITION_PLACEHOLDER
#SBATCH --output=LOG_DIR_PLACEHOLDER/slurm_batch_BATCH_ID_PLACEHOLDER_%j.out
#SBATCH --error=LOG_DIR_PLACEHOLDER/slurm_batch_BATCH_ID_PLACEHOLDER_%j.err

set -e

echo "Job started at: $(date)"
echo "Processing protein: PROTEIN_PLACEHOLDER, batch: BATCH_ID_PLACEHOLDER"

# Activate environment
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate RFAA_CONDA_ENV_PLACEHOLDER

# Set database paths for MSA generation
export DB_UR30="DB_UR30_PLACEHOLDER"
export DB_BFD="DB_BFD_PLACEHOLDER"

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
BATCH_TARGETS_FILE="${LOCAL_OUT}/batch_targets.txt"
: > "$BATCH_TARGETS_FILE"

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
    echo "compound_${COMPOUND_ID}" >> "$BATCH_TARGETS_FILE"

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

database_params:
  hhdb: /home/yangl_pacagen_com/Applications/model_weights/rosetta_db/pdb100_2021Mar03/pdb100_2021Mar03

protein_inputs:
  A:
    fasta_file: "FASTA_FILE_PLACEHOLDER"

sm_inputs:
  B:
    input: |-
      ${SMILES}
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

# Copy results directly to output directory
if [ $SUCCESS_COUNT -gt 0 ]; then
    echo "Copying results..."

    for compound_dir in "$LOCAL_OUT"/compound_*; do
        if [ ! -d "$compound_dir" ]; then
            continue
        fi

        compound_name=$(basename "$compound_dir")
        compound_id="${compound_name#compound_}"
        if [ "$compound_name" = "$compound_id" ] || [ -z "$compound_id" ]; then
            echo "Warning: Unexpected compound directory name: $compound_dir"
            continue
        fi

        output_dir="OUTPUT_DIR_PLACEHOLDER/compound_${compound_id}"
        mkdir -p "$output_dir"

        # Extract and save SMILES from config file (but don't copy the config)
        config_file="$LOCAL_OUT/config_${compound_name}.yaml"
        if [ -f "$config_file" ]; then
            python3 -c "
import yaml
with open('$config_file') as f:
    cfg = yaml.safe_load(f)
smiles = cfg.get('sm_inputs', {}).get('B', {}).get('input', '')
with open('$output_dir/smiles.txt', 'w') as f:
    f.write(smiles)
" 2>/dev/null || echo "Warning: Could not extract SMILES for compound $compound_id"
        fi

        # Copy aux.pt file (contains scores: mean_plddt, mean_pae, pae_prot, pae_inter)
        aux_file=$(find "$compound_dir" -maxdepth 1 -name "*_aux.pt" | head -1)
        if [ -f "$aux_file" ]; then
            cp "$aux_file" "$output_dir/"
            echo "  Copied aux.pt for compound $compound_id"
        fi

        # Copy best PDB structure (typically only one main model)
        pdb_file=$(find "$compound_dir" -maxdepth 1 -name "*.pdb" | head -1)
        if [ -f "$pdb_file" ]; then
            cp "$pdb_file" "$output_dir/structure.pdb"
        fi
    done

    echo "Results copied to OUTPUT_DIR_PLACEHOLDER/"
fi

# Validate copied outputs and append batch result to job ledger
EXPECTED_COUNT=$(sort -u "$BATCH_TARGETS_FILE" | sed '/^$/d' | wc -l)
COMPLETED_COUNT=0
while IFS= read -r target_compound; do
    [ -z "$target_compound" ] && continue
    if compgen -G "OUTPUT_DIR_PLACEHOLDER/${target_compound}/*_aux.pt" > /dev/null; then
        COMPLETED_COUNT=$((COMPLETED_COUNT + 1))
    fi
done < <(sort -u "$BATCH_TARGETS_FILE")

JOB_STATUS="failed"
if [ "$EXPECTED_COUNT" -gt 0 ] && [ "$COMPLETED_COUNT" -eq "$EXPECTED_COUNT" ]; then
    JOB_STATUS="success"
fi

source "LEDGER_HELPER_PLACEHOLDER"
append_job_result "LEDGER_FILE_PLACEHOLDER" "${SLURM_JOB_ID:-unknown}" "PROTEIN_PLACEHOLDER" "batch_BATCH_ID_PLACEHOLDER" "$JOB_STATUS" "$COMPLETED_COUNT" "$EXPECTED_COUNT"
echo "Batch ledger status: ${JOB_STATUS} (${COMPLETED_COUNT}/${EXPECTED_COUNT})"

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
    sed -i "s|PARTITION_PLACEHOLDER|${LIGAND_PARTITION}|g" "$JOB_SCRIPT"
    sed -i "s|LOG_DIR_PLACEHOLDER|${LOG_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|RFAA_CONDA_ENV_PLACEHOLDER|${RFAA_CONDA_ENV}|g" "$JOB_SCRIPT"
    sed -i "s|RFAA_ROOT_PLACEHOLDER|${RFAA_ROOT}|g" "$JOB_SCRIPT"
    sed -i "s|SELECTED_CSV_PLACEHOLDER|${SELECTED_CSV}|g" "$JOB_SCRIPT"
    sed -i "s|N_BATCHES_PLACEHOLDER|${N_BATCHES}|g" "$JOB_SCRIPT"
    sed -i "s|FASTA_FILE_PLACEHOLDER|${FASTA_FILE}|g" "$JOB_SCRIPT"
    sed -i "s|OUTPUT_DIR_PLACEHOLDER|${OUTPUT_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|TOKEN_FILE_PLACEHOLDER|${TOKEN_FILE}|g" "$JOB_SCRIPT"
    sed -i "s|LEDGER_FILE_PLACEHOLDER|${LEDGER_FILE}|g" "$JOB_SCRIPT"
    sed -i "s|LEDGER_HELPER_PLACEHOLDER|${LEDGER_HELPER}|g" "$JOB_SCRIPT"
    sed -i "s|PROTEIN_FOLD_OUTPUT_PLACEHOLDER|${PROTEIN_FOLD_OUTPUT}|g" "$JOB_SCRIPT"
    sed -i "s|DB_UR30_PLACEHOLDER|${ROSETTA_DB_UR30}|g" "$JOB_SCRIPT"
    sed -i "s|DB_BFD_PLACEHOLDER|${ROSETTA_DB_BFD}|g" "$JOB_SCRIPT"

    # Submit job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted protein-ligand batch ${batch_id} for ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${RFAA_DIR}/protein_ligand/job_ids.txt"
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
    mkdir -p "${RFAA_DIR}/protein_ligand"
    mkdir -p "${RFAA_DIR}/protein_ligand/output"
    > "${RFAA_DIR}/protein_ligand/job_ids.txt"
    > "${RFAA_DIR}/protein_ligand/output/job_results.tsv"

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
