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
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

# Load configuration from environment or use defaults
TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
SCRIPT_ROOT="/home/yangl_pacagen_com/screening_workflow/scripts"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"

# Tools
VINA_EXE="/home/yangl_pacagen_com/screening_workflow/scripts/docking/docking.py"

# SLURM configuration
TIME_LIMIT="48:00:00"     # 8 hours per job
MEMORY="15G"              # Memory per job
CPUS_PER_TASK=2

# GPU configuration (Vina might not need GPU, adjust as needed)
# Set to empty string if no GPU needed
GPU_REQUEST=""            # Empty = no GPU request
# GPU_REQUEST="--gres=gpu:1"  # Uncomment if GPU is needed

# SLURM partition
PARTITION="g24"

# Receptor prep job id (set by submit_receptor_prep_job)
RECEPTOR_PREP_JOB_ID=""

# Poll interval while waiting for receptor prep job completion
PREP_POLL_INTERVAL=30

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
import json

df = pd.read_csv("${seq_file}")
row = df.loc[df['protein_name'] == "${protein}"].iloc[0]

cx = float(row['center_x'])
cy = float(row['center_y'])
cz = float(row['center_z'])
sx = float(row['size_x'])
sy = float(row['size_y'])
sz = float(row['size_z'])

print(json.dumps([[cx, cy, cz, sx, sy, sz]]))
EOF
}

is_receptor_cache_ready() {
    local receptor_dir=$1
    local boxes_json=$2

    python3 - "$receptor_dir" "$boxes_json" <<'PY'
import json
import os
import sys

receptor_dir = sys.argv[1]
boxes_json = sys.argv[2]
manifest_path = os.path.join(receptor_dir, "receptor_manifest.json")

if not os.path.exists(manifest_path):
    sys.exit(1)

def normalize_boxes(values):
    if values and not isinstance(values[0], (list, tuple)):
        values = [values]
    return [[float(x) for x in row] for row in values]

def close_boxes(a, b, tol=1e-3):
    if len(a) != len(b):
        return False
    for row_a, row_b in zip(a, b):
        if len(row_a) != len(row_b):
            return False
        for x, y in zip(row_a, row_b):
            if abs(float(x) - float(y)) > tol:
                return False
    return True

try:
    with open(manifest_path, "r") as f:
        manifest = json.load(f)

    manifest_boxes = normalize_boxes(manifest["boxes_input"])
    requested_boxes = normalize_boxes(json.loads(boxes_json))

    if not close_boxes(manifest_boxes, requested_boxes):
        sys.exit(1)

    receptor_entries = manifest["receptor_pdbqts"]
    if isinstance(receptor_entries, dict):
        receptor_entries = [receptor_entries[k] for k in sorted(receptor_entries, key=lambda x: int(x))]

    if len(receptor_entries) != len(manifest_boxes):
        sys.exit(1)

    for entry in receptor_entries:
        receptor_path = entry
        if not os.path.isabs(receptor_path):
            receptor_path = os.path.join(receptor_dir, receptor_path)
        if not os.path.exists(receptor_path) or os.path.getsize(receptor_path) == 0:
            sys.exit(1)
except Exception:
    sys.exit(1)

sys.exit(0)
PY
}

submit_receptor_prep_job() {
    local protein=$1
    local boxes_json=$2

    local FINE_DIR="${TASK_ROOT}/${protein}/fine_screening"
    local OUTPUT_DIR="${FINE_DIR}/Vina/output"
    local RECEPTOR_DIR="${FINE_DIR}/Vina/receptor"
    local PDB_FILE="${TASK_ROOT}/Input/protein_file/${protein}/${protein}.pdb"
    local JOB_SCRIPT="${OUTPUT_DIR}/slurm_receptor_prep.sh"

    RECEPTOR_PREP_JOB_ID=""

    if [ ! -f "$PDB_FILE" ]; then
        log_error "PDB file not found: $PDB_FILE"
        return 1
    fi

    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$RECEPTOR_DIR"

    if is_receptor_cache_ready "$RECEPTOR_DIR" "$boxes_json"; then
        log_info "Reusing prepared receptor cache for ${protein}: ${RECEPTOR_DIR}"
        return 0
    fi

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=vina_prep_${protein}
EOF

    # Add GPU request if specified
    if [ -n "$GPU_REQUEST" ]; then
        echo "#SBATCH ${GPU_REQUEST}" >> "$JOB_SCRIPT"
    fi

    cat >> "$JOB_SCRIPT" <<EOF
#SBATCH --time=${TIME_LIMIT}
#SBATCH --mem=${MEMORY}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --partition=${PARTITION}
#SBATCH --output=${OUTPUT_DIR}/slurm_receptor_prep_%j.out
#SBATCH --error=${OUTPUT_DIR}/slurm_receptor_prep_%j.err

set -e

echo "Receptor prep job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Protein: ${protein}"
echo "Docking boxes: ${boxes_json}"

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate vina_test

mkdir -p "${RECEPTOR_DIR}"

python "${VINA_EXE}" \\
    --pdb "${PDB_FILE}" \\
    --boxes '${boxes_json}' \\
    --output "${RECEPTOR_DIR}" \\
    --prepare-receptor-only

echo "Receptor prep completed at: \$(date)"
EOF

    local job_id
    job_id=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$job_id" ]; then
        RECEPTOR_PREP_JOB_ID="$job_id"
        log_info "Submitted receptor prep for ${protein} (Job ID: ${job_id})"
        echo "$job_id" >> "${OUTPUT_DIR}/receptor_job_ids.txt"
        return 0
    fi

    log_error "Failed to submit receptor prep for ${protein}"
    return 1
}

wait_for_receptor_prep_job() {
    local protein=$1
    local job_id=$2
    local receptor_dir=$3
    local boxes_json=$4

    if [ -z "$job_id" ]; then
        return 0
    fi

    log_info "Waiting for receptor prep job ${job_id} (${protein}) to complete..."

    while squeue -j "$job_id" -h >/dev/null 2>&1 && [ -n "$(squeue -j "$job_id" -h 2>/dev/null)" ]; do
        sleep "$PREP_POLL_INTERVAL"
    done

    if is_receptor_cache_ready "$receptor_dir" "$boxes_json"; then
        log_info "Receptor prep completed successfully for ${protein}"
        return 0
    fi

    log_error "Receptor prep failed for ${protein} (job ${job_id})."
    return 1
}

# Function to submit a single Vina docking job
submit_vina_job() {
    local protein=$1
    local part=$2
    local boxes_json=$3
    local receptor_dir=$4
    local prep_job_id=$5

    local FINE_DIR="${TASK_ROOT}/${protein}/fine_screening"
    local INPUT_CSV="${FINE_DIR}/Vina/input/${part}.csv"
    local OUTPUT_DIR="${FINE_DIR}/Vina/output/${part}"
    local TOKEN_FILE="${FINE_DIR}/Vina/output/${part}.done"
    local LEDGER_FILE="${FINE_DIR}/Vina/output/job_results.tsv"

    # Check if input CSV exists
    if [ ! -f "$INPUT_CSV" ]; then
        log_error "Input CSV not found: $INPUT_CSV"
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

    # Add GPU request if specified
    if [ -n "$GPU_REQUEST" ]; then
        echo "#SBATCH ${GPU_REQUEST}" >> "$JOB_SCRIPT"
    fi

    cat >> "$JOB_SCRIPT" <<EOF
#SBATCH --time=${TIME_LIMIT}
#SBATCH --mem=${MEMORY}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --partition=${PARTITION}
#SBATCH --output=${FINE_DIR}/Vina/output/slurm_${part}_%j.out
#SBATCH --error=${FINE_DIR}/Vina/output/slurm_${part}_%j.err

# Error handling
set -e

# Log start time
echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Processing part ${part} for protein ${protein}"
echo "Docking boxes: ${boxes_json}"
echo "Prepared receptor directory: ${receptor_dir}"

# Activate conda environment
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate vina_test
source "${SCRIPT_ROOT}/pipeline/lib/job_result_ledger.sh"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Run Vina docking
echo "Starting Vina docking..."
VINA_EXIT_CODE=0
python "${VINA_EXE}" \\
    --smiles "${INPUT_CSV}" \\
    --boxes '${boxes_json}' \\
    --output "${OUTPUT_DIR}" \\
    --prepared-receptor-dir "${receptor_dir}" \\
    --smiles-col "ligand_description" \\
    --skip-docked || VINA_EXIT_CODE=\$?

# Validate outputs and append part result to job ledger
EXPECTED_COUNT=\$(python3 -c "import csv; import sys; f='${INPUT_CSV}'; rows=max(sum(1 for _ in open(f, 'r', encoding='utf-8'))-1,0); print(rows)")
COMPLETED_COUNT=\$(find "${OUTPUT_DIR}" -type f -name "docking_affinities.txt" 2>/dev/null | wc -l)

JOB_STATUS="failed"
if [ "\$VINA_EXIT_CODE" -eq 0 ] && [ "\$EXPECTED_COUNT" -gt 0 ] && [ "\$COMPLETED_COUNT" -ge "\$EXPECTED_COUNT" ]; then
    JOB_STATUS="success"
fi

append_job_result "${LEDGER_FILE}" "\${SLURM_JOB_ID:-unknown}" "${protein}" "${part}" "\$JOB_STATUS" "\$COMPLETED_COUNT" "\$EXPECTED_COUNT"
echo "Part ledger status: \${JOB_STATUS} (\${COMPLETED_COUNT}/\${EXPECTED_COUNT})"

if [ "\$JOB_STATUS" = "success" ]; then
    echo "Vina docking completed successfully"
    touch "${TOKEN_FILE}"
else
    echo "Vina docking failed"
    exit 1
fi

echo "Job completed at: \$(date)"
EOF

    local sbatch_cmd=(sbatch --parsable)
    if [ -n "$prep_job_id" ]; then
        sbatch_cmd+=(--dependency="afterok:${prep_job_id}")
    fi

    # Submit the job
    local job_id
    job_id=$("${sbatch_cmd[@]}" "$JOB_SCRIPT")

    if [ -n "$job_id" ]; then
        log_info "Submitted part ${part} for protein ${protein} (Job ID: ${job_id})"
        echo "$job_id" >> "${FINE_DIR}/Vina/output/job_ids.txt"
        return 0
    fi

    log_error "Failed to submit part ${part} for protein ${protein}"
    return 1
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

    # Get docking box parameters as JSON
    BOXES_JSON=$(get_box_params "$PROTEIN" "$SEQS_CSV")
    if [ -z "$BOXES_JSON" ]; then
        log_error "Failed to extract docking box parameters for ${PROTEIN}"
        continue
    fi

    # Clear previous job lists
    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/output"
    RECEPTOR_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/receptor"
    mkdir -p "$OUTPUT_DIR"
    > "${OUTPUT_DIR}/job_ids.txt"
    > "${OUTPUT_DIR}/receptor_job_ids.txt"
    > "${OUTPUT_DIR}/job_results.tsv"

    if ! submit_receptor_prep_job "$PROTEIN" "$BOXES_JSON"; then
        log_error "Skipping ${PROTEIN}: receptor prep submission failed"
        continue
    fi

    PREP_JOB_ID="$RECEPTOR_PREP_JOB_ID"
    if [ -n "$PREP_JOB_ID" ]; then
        if ! wait_for_receptor_prep_job "$PROTEIN" "$PREP_JOB_ID" "$RECEPTOR_DIR" "$BOXES_JSON"; then
            log_error "Stopping Vina submissions because receptor prep failed for ${PROTEIN}"
            exit 1
        fi
    fi

    # Submit jobs for each part only after receptor prep succeeds
    SUBMITTED=0
    for INPUT_FILE in "${INPUT_FILES[@]}"; do
        # Extract part name (e.g., "part_0" from "/path/to/part_0.csv")
        PART=$(basename "$INPUT_FILE" .csv)

        if submit_vina_job "$PROTEIN" "$PART" "$BOXES_JSON" "$RECEPTOR_DIR" "$PREP_JOB_ID"; then
            SUBMITTED=$((SUBMITTED + 1))
        fi
    done

    if [ -n "$PREP_JOB_ID" ]; then
        log_info "Receptor prep job ID for ${PROTEIN}: ${PREP_JOB_ID}"
    fi
    log_info "Submitted ${SUBMITTED} docking jobs for protein ${PROTEIN}"
    log_info "Docking job IDs saved to: ${OUTPUT_DIR}/submitted_jobs.txt"
    log_info "Receptor prep job IDs saved to: ${OUTPUT_DIR}/receptor_job_ids.txt"
done

log_info "All docking jobs submitted!"
log_info ""
log_info "Monitor jobs with: squeue -u \$USER"
log_info "Check job outputs in: \${TASK_ROOT}/\${PROTEIN}/fine_screening/Vina/output/"
