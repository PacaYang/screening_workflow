#!/bin/bash
# Step 4: Submit AF3 SLURM jobs for peptide validation

set -euo pipefail

# Parse arguments
OUTPUT_DIR=""
PROTEIN_NAME=""
N_BATCHES=40

while [[ $# -gt 0 ]]; do
    case $1 in
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --protein-name) PROTEIN_NAME="$2"; shift 2 ;;
        --n-batches) N_BATCHES="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$OUTPUT_DIR" || -z "$PROTEIN_NAME" ]]; then
    echo "Error: Missing required arguments"
    exit 1
fi

# Setup paths
INPUT_DIR="${OUTPUT_DIR}/03_af3_input"
AF3_DIR="${OUTPUT_DIR}/04_af3_output"
OUTPUT_AF3="${AF3_DIR}/predictions"
LOG_DIR="${AF3_DIR}/logs"
TOKEN_DIR="${OUTPUT_AF3}/token"

mkdir -p "$OUTPUT_AF3" "$LOG_DIR" "$TOKEN_DIR"

# AF3 configuration
AF3_EXE="/home/yangl_pacagen_com/Applications/alphafold3/run_alphafold.py"
AF3_WEIGHT_DIR="/home/yangl_pacagen_com/Applications/model_weights/AF3"
AF3_DB_DIR="/home/yangl_pacagen_com/Applications/model_weights/af3_db"

# SLURM configuration
TIME_LIMIT="48:00:00"
MEMORY="15G"
CPUS_PER_TASK=2
PARTITION="g24"

# Check input files
N_JSONS=$(find "$INPUT_DIR" -name "*.json" | wc -l)
if [[ $N_JSONS -eq 0 ]]; then
    echo "Error: No JSON files found in $INPUT_DIR"
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found ${N_JSONS} JSON files to process"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Splitting into ${N_BATCHES} batches"

# Split JSONs into batches
JSON_LIST=($(find "$INPUT_DIR" -name "*.json" | sort))
BATCH_SIZE=$(( (N_JSONS + N_BATCHES - 1) / N_BATCHES ))

# Submit batch jobs
BATCH_ID=0
for ((i=0; i<N_JSONS; i+=BATCH_SIZE)); do
    BATCH_ID=$((BATCH_ID + 1))
    BATCH_JSONS=("${JSON_LIST[@]:i:BATCH_SIZE}")

    # Create batch input list
    BATCH_FILE="${AF3_DIR}/batch_${BATCH_ID}.txt"
    printf "%s\n" "${BATCH_JSONS[@]}" > "$BATCH_FILE"

    # Submit SLURM job
    JOB_NAME="${PROTEIN_NAME}_af3_batch${BATCH_ID}"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${LOG_DIR}/batch_${BATCH_ID}_%j.out
#SBATCH --error=${LOG_DIR}/batch_${BATCH_ID}_%j.err
#SBATCH --time=${TIME_LIMIT}
#SBATCH --mem=${MEMORY}
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --partition=${PARTITION}

set -e

echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Processing batch ${BATCH_ID}"

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate af3_test

export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

BATCH_LOCAL_OUT=/tmp/af3_batch_${PROTEIN_NAME}_${BATCH_ID}_\${SLURM_JOB_ID}
mkdir -p \$BATCH_LOCAL_OUT

SUCCESS_COUNT=0
FAIL_COUNT=0

while IFS= read -r json_file; do
    json_name=\$(basename "\$json_file" .json)
    token_file="${TOKEN_DIR}/\${json_name}.token"

    if [[ -f "\$token_file" ]]; then
        echo "Skipping \$json_name (already processed)"
        continue
    fi

    echo "Processing \$json_name"

    LOCAL_OUT=\$BATCH_LOCAL_OUT/\${json_name}
    mkdir -p \$LOCAL_OUT

    if python ${AF3_EXE} \\
        --json_path="\$json_file" \\
        --model_dir="${AF3_WEIGHT_DIR}" \\
        --db_dir="${AF3_DB_DIR}" \\
        --output_dir=\$LOCAL_OUT \\
        --norun_data_pipeline; then

        # Copy results to shared storage
        if [ -d "\$LOCAL_OUT/\${json_name}" ]; then
            cp -r "\$LOCAL_OUT/\${json_name}" "${OUTPUT_AF3}/"
            touch "\$token_file"
            SUCCESS_COUNT=\$((SUCCESS_COUNT + 1))
            echo "  ✓ \$json_name completed"
        else
            FAIL_COUNT=\$((FAIL_COUNT + 1))
            echo "  ✗ \$json_name failed (no output)"
        fi
    else
        FAIL_COUNT=\$((FAIL_COUNT + 1))
        echo "  ✗ \$json_name failed"
    fi

    rm -rf "\$LOCAL_OUT/\${json_name}"
done < "${BATCH_FILE}"

rm -rf \$BATCH_LOCAL_OUT

echo ""
echo "Batch ${BATCH_ID} summary:"
echo "  Successful: \$SUCCESS_COUNT"
echo "  Failed: \$FAIL_COUNT"
echo "Job completed at: \$(date)"
EOF

done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitted ${BATCH_ID} batch jobs"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Monitor with: squeue -u \$USER"
