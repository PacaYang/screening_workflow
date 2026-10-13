#!/bin/bash
#
# Standalone script for Snakemake rule: run_graphdta
# Submits SLURM jobs for GraphDTA initial screening (one per input chunk).
#

set -e

# ============================================================================
# Configuration
# ============================================================================

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
EXE="/home/yangl_pacagen_com/screening_workflow/scripts/initial_screening/GraphDTA_predict.py"
MODEL="${MASTER_GRAPHDTA_MODEL:-/home/yangl_pacagen_com/Applications/model_weights/GraphDTA/model_GINConvNet_kiba.pt}"
CONDA_ENV="graphdta_test"

# SLURM configuration
TIME_LIMIT="48:00:00"
MEMORY="15G"
CPUS_PER_TASK=2
PARTITION="g24"
MAX_CONCURRENT=40

# ============================================================================
# Functions
# ============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

get_proteins() {
    if [ -n "$MASTER_PROTEINS" ]; then
        echo "$MASTER_PROTEINS"
    else
        echo "JAK1JH1"
    fi
}

# Submit all pending input chunks for a protein as a single SLURM job array.
# Each array task maps its SLURM_ARRAY_TASK_ID (1-based) to an input index via
# a manifest file, so non-sequential indices are supported.
submit_array() {
    local protein=$1
    local OUTPUT_DIR=$2
    shift 2
    local indices=("$@")

    local MANIFEST="${OUTPUT_DIR}/array_manifest.txt"
    printf '%s\n' "${indices[@]}" > "$MANIFEST"
    local n=${#indices[@]}

    local JOB_SCRIPT="${OUTPUT_DIR}/slurm_array.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=graphdta_${protein}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --mem=${MEMORY}
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --partition=${PARTITION}
#SBATCH --array=1-${n}%${MAX_CONCURRENT}
#SBATCH --output=${OUTPUT_DIR}/slurm_%A_%a.out
#SBATCH --error=${OUTPUT_DIR}/slurm_%A_%a.err

set -e

echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Array task: \${SLURM_ARRAY_TASK_ID}"

IDX=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")
if [ -z "\$IDX" ]; then
    echo "ERROR: no manifest entry for array task \${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

INPUT_CSV="${TASK_ROOT}/${protein}/initial_screening/inputs/input_\${IDX}.csv"
OUTPUT_CSV="${OUTPUT_DIR}/prediction_\${IDX}.csv"
FAILED_CSV="${OUTPUT_DIR}/failed_smiles_\${IDX}.csv"

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

mkdir -p "${OUTPUT_DIR}"

export PYTHONPATH="/home/yangl_pacagen_com/Applications/GraphDTA:\${PYTHONPATH}"
python "${EXE}" \\
    "\${INPUT_CSV}" \\
    --model "${MODEL}" \\
    --model-type "gin" \\
    --batch-size "64" \\
    --output "\${OUTPUT_CSV}" \\
    --failed-smiles "\${FAILED_CSV}" \\
    --smiles-column "SMILES"

echo "Job completed at: \$(date)"
EOF

    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted array of ${n} tasks for ${protein} (Job ID: ${JOB_ID}, max ${MAX_CONCURRENT} concurrent)"
        echo "$JOB_ID" >> "${OUTPUT_DIR}/job_ids.txt"
        return 0
    else
        log_error "Failed to submit array for ${protein}"
        return 1
    fi
}

# ============================================================================
# Main
# ============================================================================

log_info "Starting GraphDTA batch job submission"

PROTEINS=$(get_proteins)

for PROTEIN in $PROTEINS; do
    log_info "Processing protein: $PROTEIN"

    INPUT_DIR="${TASK_ROOT}/${PROTEIN}/initial_screening/inputs"
    TOKEN="${INPUT_DIR}/finish.token"

    if [ ! -f "$TOKEN" ]; then
        log_error "finish.token not found for ${PROTEIN}: ${TOKEN}"
        log_error "Please run run_make_input.sh first"
        continue
    fi

    INPUT_FILES=$(ls "${INPUT_DIR}"/input_*.csv 2>/dev/null | sort)

    if [ -z "$INPUT_FILES" ]; then
        log_error "No input CSV files found in ${INPUT_DIR}"
        continue
    fi

    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/initial_screening/GraphDTA"
    mkdir -p "$OUTPUT_DIR"
    > "${OUTPUT_DIR}/job_ids.txt"

    # Collect indices whose predictions don't already exist
    PENDING=()
    for INPUT_FILE in $INPUT_FILES; do
        IDX=$(basename "$INPUT_FILE" .csv | sed 's/input_//')
        if [ -f "${OUTPUT_DIR}/prediction_${IDX}.csv" ]; then
            log_info "prediction_${IDX}.csv already exists for ${PROTEIN}, skipping"
            continue
        fi
        PENDING+=("$IDX")
    done

    if [ ${#PENDING[@]} -eq 0 ]; then
        log_info "All predictions already exist for ${PROTEIN}, nothing to submit"
        continue
    fi

    submit_array "$PROTEIN" "$OUTPUT_DIR" "${PENDING[@]}"

    log_info "Submitted ${#PENDING[@]} GraphDTA tasks for ${PROTEIN}"
done

log_info "All GraphDTA jobs submitted"
log_info "Monitor jobs with: squeue -u \$USER"
