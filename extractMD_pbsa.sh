#!/bin/bash
#
# Batch script to extract MD results and run PBSA using SLURM
# This script divides the work into 20 batches and submits them to SLURM
#
# Usage: ./extract_md_and_run_pbsa_batch.sh [base_dir]
#
# Example: ./extract_md_and_run_pbsa_batch.sh /shared/B3/IL4_IL13_IL13RA1/IL4/fine_screening/PBSA/PBSA
#

set -e

# ============================================================================
# Configuration
# ============================================================================

# Base directory containing MD subdirectories
BASE_DIR="${1:-/shared/B3/IL4_IL13_IL13RA1/IL4/fine_screening/PBSA/PBSA}"
MD_DIR="${BASE_DIR}/MD"
PBSA_DIR="${BASE_DIR}/PBSA"
LOG_DIR="${BASE_DIR}/logs"

# Number of batches to divide work into
N_BATCHES=20

# Time point to extract (in ps)
EXTRACT_TIME=4500

# PBSA configuration
PBSA_EXE="/home/ubuntu/miniconda3/envs/gmxMMPBSA/bin/gmx_MMPBSA"
GMX_RC="/home/ubuntu/Applications/gromacs-2025.3/bin/GMXRC"
PBSA_SCRIPT_DIR="/home/ubuntu/screening_workflow/scripts/pbsa"

# SLURM configuration
TIME_LIMIT="48:00:00"     # 48 hours per batch
CPUS_PER_TASK=32

# EC2 instance constraint (if using AWS ParallelCluster)
CONSTRAINT="g5.16xlarge"    # Set to empty string if not using constraints
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

# Function to submit a batch job
submit_batch_job() {
    local batch_id=$1
    local start_idx=$2
    local end_idx=$3
    shift 3
    local compound_list=("$@")

    local JOB_SCRIPT="${LOG_DIR}/slurm_batch_${batch_id}.sh"

    # Write job header
    cat > "$JOB_SCRIPT" <<EOFHEADER
#!/bin/bash
#SBATCH --job-name=md_pbsa_extract_b${batch_id}
EOFHEADER

    # Add constraint if specified
    if [ -n "$CONSTRAINT" ]; then
        echo "#SBATCH --constraint=${CONSTRAINT}" >> "$JOB_SCRIPT"
    fi

    cat >> "$JOB_SCRIPT" <<EOFHEADER
#SBATCH --time=${TIME_LIMIT}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --output=${LOG_DIR}/slurm_batch_${batch_id}_%j.out
#SBATCH --error=${LOG_DIR}/slurm_batch_${batch_id}_%j.err

# Log start time
echo "=========================================="
echo "MD Extraction + PBSA Batch Job"
echo "=========================================="
echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Processing batch ${batch_id}"
echo ""

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate gmxMMPBSA

EOFHEADER

    cat >> "$JOB_SCRIPT" <<'EOFCONFIG'

# ============================================================================
# Configuration
# ============================================================================

MD_DIR="PLACEHOLDER_MD_DIR"
PBSA_DIR="PLACEHOLDER_PBSA_DIR"
EXTRACT_TIME=PLACEHOLDER_EXTRACT_TIME
PBSA_EXE="PLACEHOLDER_PBSA_EXE"
GMX_RC="PLACEHOLDER_GMX_RC"
PBSA_SCRIPT_DIR="PLACEHOLDER_PBSA_SCRIPT_DIR"

# ============================================================================
# Functions
# ============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}
EOFCONFIG

    # Replace configuration placeholders
    sed -i "s|PLACEHOLDER_MD_DIR|${MD_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|PLACEHOLDER_PBSA_DIR|${PBSA_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|PLACEHOLDER_EXTRACT_TIME|${EXTRACT_TIME}|g" "$JOB_SCRIPT"
    sed -i "s|PLACEHOLDER_PBSA_EXE|${PBSA_EXE}|g" "$JOB_SCRIPT"
    sed -i "s|PLACEHOLDER_GMX_RC|${GMX_RC}|g" "$JOB_SCRIPT"
    sed -i "s|PLACEHOLDER_PBSA_SCRIPT_DIR|${PBSA_SCRIPT_DIR}|g" "$JOB_SCRIPT"

    cat >> "$JOB_SCRIPT" <<'EOFFUNCTIONS'

# Function to extract MD frame and run PBSA for a single compound
process_compound() {
    local compound_name=$1

    log_info "Processing ${compound_name}"

    local compound_dir="${MD_DIR}/${compound_name}"

    # Paths for this compound
    local xtc_file="${compound_dir}/T298.xtc"
    local tpr_file="${compound_dir}/T298.tpr"
    local gro_file="${compound_dir}/T298.gro"
    local top_file="${compound_dir}/system.top"
    local ndx_file="${compound_dir}/index.ndx"

    local pbsa_outdir="${PBSA_DIR}/${compound_name}"
    local pbsa_token="${pbsa_outdir}/token.done"

    # Check if PBSA already completed
    if [ -f "$pbsa_token" ]; then
        log_info "  SKIPPING: PBSA already completed for ${compound_name}"
        return 0
    fi

    # Check if compound directory exists
    if [ ! -d "$compound_dir" ]; then
        log_error "  SKIPPING: Compound directory not found: $compound_dir"
        return 1
    fi

    # Check if required MD files exist
    if [ ! -f "$xtc_file" ]; then
        log_error "  SKIPPING: XTC file not found: $xtc_file"
        return 1
    fi

    if [ ! -f "$tpr_file" ]; then
        log_error "  SKIPPING: TPR file not found: $tpr_file"
        return 1
    fi

    if [ ! -f "$top_file" ]; then
        log_error "  SKIPPING: TOP file not found: $top_file"
        return 1
    fi

    if [ ! -f "$ndx_file" ]; then
        log_error "  SKIPPING: NDX file not found: $ndx_file"
        return 1
    fi
    
    # ========================================================================
    # STEP 2: Run PBSA Analysis
    # ========================================================================
    log_info "  Running PBSA analysis"

    # Set up environment for PBSA
    export PATH="/home/ubuntu/miniconda3/envs/gmxMMPBSA/bin:$PATH"

    # Create PBSA output directory
    mkdir -p "$pbsa_outdir"
    cd "$pbsa_outdir"

    # Output files
    local dat_file="${pbsa_outdir}/FINAL_RESULTS_MMPBSA.dat"
    local csv_file="${pbsa_outdir}/FINAL_RESULTS_MMPBSA.csv"

    # Run PBSA calculation
    if "$PBSA_EXE" -O \
        -i "${PBSA_SCRIPT_DIR}/mmpbsa.in" \
        -cs "$tpr_file" \
        -ct "$xtc_file" \
        -ci "$ndx_file" \
        -cg 1 13 \
        -cp "$top_file" \
        -o "$dat_file" \
        -eo "$csv_file"; then

        log_info "  PBSA completed successfully"
        touch "$pbsa_token"
        return 0
    else
        log_error "  PBSA failed"
        return 1
    fi
}
EOFFUNCTIONS

    # Write compound list header
    cat >> "$JOB_SCRIPT" <<'EOFLIST'

# ============================================================================
# Main Processing Loop
# ============================================================================

# List of compounds to process in this batch
COMPOUNDS=(
EOFLIST

    # Write compound list directly to the file
    for ((i=$start_idx; i<$end_idx; i++)); do
        echo "    \"${compound_list[$i]}\"" >> "$JOB_SCRIPT"
    done

    # Continue writing the rest of the script
    cat >> "$JOB_SCRIPT" <<'EOFMAIN'
)

echo "Batch processing: ${#COMPOUNDS[@]} compounds"
echo ""

SUCCESS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

for compound_name in "${COMPOUNDS[@]}"; do
    echo "=========================================="
    if process_compound "$compound_name"; then
        if [ -f "${PBSA_DIR}/${compound_name}/token.done" ]; then
            if [ $? -eq 0 ]; then
                SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            else
                SKIP_COUNT=$((SKIP_COUNT + 1))
            fi
        else
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        fi
    else
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    echo ""
done

# Summary
echo "=========================================="
echo "Batch Summary"
echo "=========================================="
echo "Processed: ${#COMPOUNDS[@]}"
echo "Successful: $SUCCESS_COUNT"
echo "Failed: $FAIL_COUNT"
echo "Skipped: $SKIP_COUNT"
echo ""
echo "Job completed at: $(date)"
echo "=========================================="
EOFMAIN

    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted batch ${batch_id} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${BASE_DIR}/submitted_jobs.txt"
        return 0
    else
        log_error "Failed to submit batch ${batch_id}"
        return 1
    fi
}

# ============================================================================
# Main Script
# ============================================================================

log_info "Starting MD extraction and PBSA batch job submission"
log_info "Base directory: $BASE_DIR"
log_info "MD directory: $MD_DIR"
log_info "PBSA directory: $PBSA_DIR"
log_info "Number of batches: $N_BATCHES"
log_info "Extraction time: ${EXTRACT_TIME} ps"
log_info ""

# Check if MD directory exists
if [ ! -d "$MD_DIR" ]; then
    log_error "MD directory not found: $MD_DIR"
    exit 1
fi

# Create output directories
mkdir -p "$LOG_DIR"
mkdir -p "$PBSA_DIR"

# Find all compound directories matching IL4_*
COMPOUND_DIRS=$(find "$MD_DIR" -maxdepth 1 -type d -name "IL4_*" | sort)

if [ -z "$COMPOUND_DIRS" ]; then
    log_error "No compound directories found matching IL4_* in $MD_DIR"
    exit 1
fi

# Build array of compound names (basename only)
COMPOUND_NAMES=()
while IFS= read -r dir; do
    COMPOUND_NAMES+=("$(basename "$dir")")
done <<< "$COMPOUND_DIRS"

N_COMPOUNDS=${#COMPOUND_NAMES[@]}
log_info "Found ${N_COMPOUNDS} compounds to process"

# Calculate batch size
BATCH_SIZE=$(( ($N_COMPOUNDS + $N_BATCHES - 1) / $N_BATCHES ))
log_info "Batch size: ${BATCH_SIZE} compounds per batch"
log_info ""

# Clear previous job list
> "${BASE_DIR}/submitted_jobs.txt"

# Submit batch jobs
SUBMITTED=0
for BATCH_ID in $(seq 0 $((N_BATCHES - 1))); do
    START_IDX=$((BATCH_ID * BATCH_SIZE))
    END_IDX=$((START_IDX + BATCH_SIZE))

    if [ $START_IDX -ge $N_COMPOUNDS ]; then
        log_info "Batch ${BATCH_ID} has no jobs, skipping"
        continue
    fi

    if [ $END_IDX -gt $N_COMPOUNDS ]; then
        END_IDX=$N_COMPOUNDS
    fi

    N_IN_BATCH=$((END_IDX - START_IDX))
    log_info "Preparing batch ${BATCH_ID}: compounds ${START_IDX} to $((END_IDX - 1)) (${N_IN_BATCH} compounds)"

    if submit_batch_job "$BATCH_ID" "$START_IDX" "$END_IDX" "${COMPOUND_NAMES[@]}"; then
        SUBMITTED=$((SUBMITTED + 1))
    fi
done

log_info ""
log_info "Submitted ${SUBMITTED} batch jobs"
log_info "Job IDs saved to: ${BASE_DIR}/submitted_jobs.txt"
log_info ""
log_info "Monitor jobs with: squeue -u \$USER"
log_info "Check job outputs in: ${LOG_DIR}/"
