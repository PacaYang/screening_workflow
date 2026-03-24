#!/bin/bash
#
# MD + PBSA Batch Automation Script
# This script submits batch jobs, each processing multiple compounds
# Grouped into 100 batches to reduce job submission overhead
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
MD_SCRIPT="${SCRIPT_DIR}/pbsa/run_pbsa_md.sh"
PBSA_EXE="/home/ubuntu/miniconda3/envs/gmxMMPBSA/bin/gmx_MMPBSA"
GMX_RC="/home/ubuntu/Applications/gromacs-2025.3/bin/GMXRC"
PBSA_SCRIPT_DIR="${SCRIPT_DIR}/pbsa"

# Number of batches to split jobs into
N_BATCHES=100

# SLURM configuration
TIME_LIMIT="240:00:00"     # 240 hours per batch
CPUS_PER_TASK=16

# EC2 instance constraint (if using AWS ParallelCluster)
CONSTRAINT="g5.4xlarge"    # Set to empty string if not using constraints
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

# Function to submit a batch MD+PBSA job
submit_md_pbsa_batch() {
    local protein=$1
    local batch_id=$2

    local FINE_DIR="${TASK_ROOT}/${protein}/fine_screening"
    local DIFFDOCK_OUTPUT="${FINE_DIR}/PBSA/DiffDock/output"
    local GRO_FILE="${TASK_ROOT}/Input/protein_file/${protein}/${protein}.gro"
    local TOP_FILE="${TASK_ROOT}/Input/protein_file/${protein}/system_EM.top"
    local LOG_DIR="${FINE_DIR}/PBSA/PBSA/logs"
    local TOKEN_DIR="${FINE_DIR}/PBSA/PBSA/token"

    # Check if required input files exist
    if [ ! -f "$GRO_FILE" ]; then
        log_error "GRO file not found: $GRO_FILE"
        return 1
    fi

    if [ ! -f "$TOP_FILE" ]; then
        log_error "TOP file not found: $TOP_FILE"
        return 1
    fi

    # Create output directories
    mkdir -p "$LOG_DIR"
    mkdir -p "$TOKEN_DIR"

    # Create SLURM job script
    local JOB_SCRIPT="${LOG_DIR}/slurm_batch_${batch_id}.sh"

    cat > "$JOB_SCRIPT" <<'EOFMAIN'
#!/bin/bash
#SBATCH --job-name=md_pbsa_PROTEIN_bBATCHID
EOFMAIN

    # Add constraint if specified
    if [ -n "$CONSTRAINT" ]; then
        echo "#SBATCH --constraint=${CONSTRAINT}" >> "$JOB_SCRIPT"
    fi

    # Add GPU request if specified
    if [ -n "$GPU_REQUEST" ]; then
        echo "#SBATCH ${GPU_REQUEST}" >> "$JOB_SCRIPT"
    fi

    cat >> "$JOB_SCRIPT" <<EOFMAIN
#SBATCH --time=${TIME_LIMIT}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --output=${LOG_DIR}/slurm_batch_${batch_id}_%j.out
#SBATCH --error=${LOG_DIR}/slurm_batch_${batch_id}_%j.err

# Log start time
echo "=========================================="
echo "MD + PBSA Batch Job"
echo "=========================================="
echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Processing batch ${batch_id} for protein ${protein}"
echo ""

# Get selected CSV to determine job range
SELECTED_CSV="${TASK_ROOT}/${protein}/initial_screening/selected.csv"
N_COMPOUNDS=\$(tail -n +2 "\$SELECTED_CSV" | wc -l)

# Calculate job range for this batch
BATCH_SIZE=\$(( (\$N_COMPOUNDS + ${N_BATCHES} - 1) / ${N_BATCHES} ))
START_IDX=\$(( ${batch_id} * \$BATCH_SIZE ))
END_IDX=\$(( \$START_IDX + \$BATCH_SIZE ))
if [ \$END_IDX -gt \$N_COMPOUNDS ]; then
    END_IDX=\$N_COMPOUNDS
fi

echo "Total compounds: \$N_COMPOUNDS"
echo "Batch size: \$BATCH_SIZE"
echo "Processing compounds \$START_IDX to \$((\$END_IDX - 1))"
echo ""

# Process each compound in the batch
SUCCESS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

for i in \$(seq \$START_IDX \$((\$END_IDX - 1))); do
    echo "=========================================="
    echo "Processing compound ${protein}_\${i} (batch ${batch_id})"
    echo "=========================================="

    SDF_FILE="${DIFFDOCK_OUTPUT}/${protein}_\${i}/rank1.sdf"
    MD_OUTDIR="${FINE_DIR}/PBSA/PBSA/MD/${protein}_\${i}"
    MD_TOKEN="\${MD_OUTDIR}/token.done"
    PBSA_OUTDIR="${FINE_DIR}/PBSA/PBSA/PBSA/${protein}_\${i}"
    PBSA_TOKEN="\${PBSA_OUTDIR}/token.done"

    # Check if SDF file exists (from DiffDock)
    if [ ! -f "\${SDF_FILE}" ]; then
        echo "SKIPPING: SDF file not found: \${SDF_FILE}"
        echo "This compound likely failed in DiffDock or was not processed."
        SKIP_COUNT=\$((SKIP_COUNT + 1))
        echo ""
        continue
    fi

    # Check if already completed
    if [ -f "\${MD_TOKEN}" ] && [ -f "\${PBSA_TOKEN}" ]; then
        echo "SKIPPING: Already completed (tokens exist)"
        SKIP_COUNT=\$((SKIP_COUNT + 1))
        echo ""
        continue
    fi

    # ========================================================================
    # STEP 1: Run Molecular Dynamics
    # ========================================================================
    echo "STEP 1: Running Molecular Dynamics"
    echo "Input SDF: \${SDF_FILE}"
    echo ""

    # Create local temporary directory on NVMe/SSD to avoid I/O contention
    LOCAL_MD_DIR=/tmp/md_${protein}_\${i}_\${SLURM_JOB_ID}
    mkdir -p \$LOCAL_MD_DIR

    echo "Using local temporary directory: \$LOCAL_MD_DIR"

    # Activate conda environment
    source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
    conda activate gmxMMPBSA

    # Run MD script in local temporary directory
    if bash "${MD_SCRIPT}" \\
        "\${SDF_FILE}" \\
        "\$LOCAL_MD_DIR" \\
        "${GRO_FILE}" \\
        "${TOP_FILE}" \\
        "${PBSA_SCRIPT_DIR}"; then

        echo "MD step completed successfully"

        # Copy results from local to network storage
        echo "Copying MD results from local storage to network storage..."
        mkdir -p "\${MD_OUTDIR}"

        if rsync -az "\$LOCAL_MD_DIR/" "\${MD_OUTDIR}/"; then
            echo "MD results copied successfully"
            touch "\${MD_TOKEN}"
        else
            echo "Failed to copy MD results"
            rm -rf "\$LOCAL_MD_DIR"
            FAIL_COUNT=\$((FAIL_COUNT + 1))
            echo ""
            continue
        fi
    else
        echo "MD step failed"
        rm -rf "\$LOCAL_MD_DIR"
        FAIL_COUNT=\$((FAIL_COUNT + 1))
        echo ""
        continue
    fi

    # Update file paths to point to network storage
    XTC_FILE="\${MD_OUTDIR}/T298.xtc"
    TPR_FILE="\${MD_OUTDIR}/T298.tpr"
    MD_TOP_FILE="\${MD_OUTDIR}/system.top"
    NDX_FILE="\${MD_OUTDIR}/index.ndx"
    DAT_FILE="\${PBSA_OUTDIR}/FINAL_RESULTS_MMPBSA.dat"

    # ========================================================================
    # STEP 2: Run PBSA Analysis
    # ========================================================================
    echo ""
    echo "STEP 2: Running PBSA Analysis"

    # Verify MD outputs exist
    if [ ! -f "\${XTC_FILE}" ] || [ ! -f "\${TPR_FILE}" ] || [ ! -f "\${MD_TOP_FILE}" ] || [ ! -f "\${NDX_FILE}" ]; then
        echo "ERROR: MD output files not found"
        rm -rf "\$LOCAL_MD_DIR"
        FAIL_COUNT=\$((FAIL_COUNT + 1))
        echo ""
        continue
    fi

    # Clean up local MD directory
    rm -rf "\$LOCAL_MD_DIR"

    # Set up environment for PBSA
    set +eu
    source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
    conda activate gmxMMPBSA
    bash -c "source ${GMX_RC}"
    export PATH="/home/ubuntu/miniconda3/envs/gmxMMPBSA/bin:\$PATH"

    # Create PBSA output directory
    mkdir -p "\${PBSA_OUTDIR}"
    cd "\${PBSA_OUTDIR}"

    # Run PBSA calculation
    if env "PATH=\$PATH" "${PBSA_EXE}" -O \\
        -i "${PBSA_SCRIPT_DIR}/mmpbsa.in" \\
        -cs "\${TPR_FILE}" \\
        -ct "\${XTC_FILE}" \\
        -ci "\${NDX_FILE}" \\
        -cg 1 13 \\
        -cp "\${MD_TOP_FILE}" \\
        -o "\${DAT_FILE}" \\
        -eo "\${PBSA_OUTDIR}/FINAL_RESULTS_MMPBSA.csv"; then

        echo "PBSA step completed successfully"
        touch "\${PBSA_TOKEN}"
        SUCCESS_COUNT=\$((SUCCESS_COUNT + 1))
    else
        echo "PBSA step failed"
        FAIL_COUNT=\$((FAIL_COUNT + 1))
    fi

    set -e
    echo ""
done

# Create completion token for this batch
touch "${TOKEN_DIR}/batch_${batch_id}.done"

echo "=========================================="
echo "Batch ${batch_id} Summary"
echo "=========================================="
echo "Successful: \$SUCCESS_COUNT"
echo "Failed: \$FAIL_COUNT"
echo "Skipped: \$SKIP_COUNT"
echo "Job completed at: \$(date)"
echo "=========================================="
EOFMAIN

    # Replace placeholders
    sed -i "s/PROTEIN/${protein}/g" "$JOB_SCRIPT"
    sed -i "s/BATCHID/${batch_id}/g" "$JOB_SCRIPT"

    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted batch ${batch_id} for protein ${protein} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${FINE_DIR}/PBSA/PBSA/submitted_jobs.txt"
        return 0
    else
        log_error "Failed to submit batch ${batch_id} for protein ${protein}"
        return 1
    fi
}

# ============================================================================
# Main Script
# ============================================================================

log_info "Starting MD + PBSA batch job submission"
log_info "Number of batches per protein: ${N_BATCHES}"

# Get list of proteins to process
PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    # Check if DiffDock output exists
    DIFFDOCK_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/DiffDock/output"
    if [ ! -d "$DIFFDOCK_DIR" ]; then
        log_error "DiffDock output not found for ${PROTEIN}: ${DIFFDOCK_DIR}"
        log_error "Please run DiffDock first"
        continue
    fi

    # Get number of compounds from selected.csv
    SELECTED_CSV="${TASK_ROOT}/${PROTEIN}/initial_screening/selected.csv"
    if [ ! -f "$SELECTED_CSV" ]; then
        log_error "Selected CSV not found: ${SELECTED_CSV}"
        continue
    fi

    N_COMPOUNDS=$(tail -n +2 "$SELECTED_CSV" | wc -l)
    BATCH_SIZE=$(( ($N_COMPOUNDS + $N_BATCHES - 1) / $N_BATCHES ))

    log_info "Found ${N_COMPOUNDS} compounds for ${PROTEIN}"
    log_info "Batch size: ${BATCH_SIZE} compounds per batch"

    # Clear previous job list
    PBSA_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/PBSA"
    mkdir -p "$PBSA_DIR"
    > "${PBSA_DIR}/submitted_jobs.txt"

    # Submit batch jobs
    SUBMITTED=0
    for BATCH_ID in $(seq 0 $((N_BATCHES - 1))); do
        # Check if this batch has jobs
        START_IDX=$((BATCH_ID * BATCH_SIZE))
        if [ $START_IDX -lt $N_COMPOUNDS ]; then
            if submit_md_pbsa_batch "$PROTEIN" "$BATCH_ID"; then
                SUBMITTED=$((SUBMITTED + 1))
            fi
        else
            log_info "Batch ${BATCH_ID} has no jobs, skipping"
        fi
    done

    log_info "Submitted ${SUBMITTED} batch jobs for protein ${PROTEIN}"
    log_info "Job IDs saved to: ${PBSA_DIR}/submitted_jobs.txt"
done

log_info "All MD + PBSA batch jobs submitted!"
log_info ""
log_info "Monitor jobs with: squeue -u \$USER"
log_info "Check job outputs in: \${TASK_ROOT}/\${PROTEIN}/fine_screening/PBSA/PBSA/logs/"
