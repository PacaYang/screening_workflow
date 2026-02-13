#!/bin/bash
#
# Control Molecules Workflow Script
# This script runs workflows (Vina, Boltz2, AlphaFold3, DiffDock, MD+PBSA) on control molecules
# Input: task_root/Input/controls.csv with columns: target, SMILES
#

set -e

# ============================================================================
# Configuration - Update these paths as needed
# ============================================================================

# Task root directory
TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"

# Script directory (where the automation scripts are located)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input files
CONTROLS_CSV="${TASK_ROOT}/Input/controls.csv"
SEQUENCES_CSV="${TASK_ROOT}/Input/sequences.csv"

# Output directory for control results
CONTROLS_OUTPUT_DIR="${TASK_ROOT}/Controls"

# Script root for helper scripts
SCRIPT_ROOT="${SCRIPT_DIR}/scripts"

# Tools and executables
VINA_EXE="${SCRIPT_ROOT}/docking.py"
BOLTZ2_EXE="/home/ubuntu/miniconda3/envs/boltz/bin/boltz"
AF3_EXE="/home/ubuntu/Applications/alphafold3/run_alphafold.py"
DIFFDOCK_DIR="/home/ubuntu/Applications/DiffDock/"
DIFFDOCK_CONFIG="/home/ubuntu/Applications/DiffDock/default_inference_args.yaml"
MD_SCRIPT="${SCRIPT_ROOT}/pbsa/run_pbsa_md.sh"
PBSA_EXE="/home/ubuntu/miniconda3/envs/gmxMMPBSA/bin/gmx_MMPBSA"
GMX_RC="/home/ubuntu/Applications/gromacs-2025.3/bin/GMXRC"
PBSA_SCRIPT_DIR="${SCRIPT_ROOT}/pbsa"

# Workflow control flags (set to 1 to enable, 0 to disable)
RUN_VINA=1
RUN_BOLTZ2=1
RUN_AF3=1
RUN_DIFFDOCK=1
RUN_MD_PBSA=1

# SLURM configuration
VINA_TIME_LIMIT="08:00:00"
VINA_MEMORY="16G"
VINA_CPUS=4

BOLTZ2_TIME_LIMIT="48:00:00"
BOLTZ2_MEMORY="15G"
BOLTZ2_CPUS=4

AF3_TIME_LIMIT="48:00:00"
AF3_MEMORY="20G"
AF3_CPUS=8

DIFFDOCK_TIME_LIMIT="48:00:00"
DIFFDOCK_MEMORY="15G"
DIFFDOCK_CPUS=4

MD_PBSA_TIME_LIMIT="240:00:00"
MD_PBSA_CPUS=64

# EC2 instance constraint
CONSTRAINT="g5.xlarge"

# ============================================================================
# Functions
# ============================================================================

log_info() {
    echo "=========================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
    echo "=========================================="
    echo ""
}

log_error() {
    echo "=========================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*"
    echo "=========================================="
    echo ""
}

# Function to extract docking box from sequences.csv
get_box_params() {
    local protein=$1

    # Use Python to parse the docking box
    python3 <<EOF
import pandas as pd
import ast

df = pd.read_csv("${SEQUENCES_CSV}")
row = df.loc[df['name'] == "${protein}"]

if row.empty:
    print("ERROR: Protein ${protein} not found in sequences.csv")
    exit(1)

cell = row['docking box'].iloc[0]

# Parse the box specification
if isinstance(cell, str):
    vals = ast.literal_eval(cell)
else:
    vals = list(cell)

center = vals[:3]
size = vals[3:]

print(" ".join(map(str, center)))
print(" ".join(map(str, size)))
EOF
}

# Function to create control-specific input directory structure
setup_control_dirs() {
    local target=$1
    local control_id=$2

    local target_dir="${CONTROLS_OUTPUT_DIR}/${target}/control_${control_id}"
    mkdir -p "${target_dir}"/{Vina,Boltz2,AF3,DiffDock,MD_PBSA}/output
    mkdir -p "${target_dir}/input"

    echo "$target_dir"
}

# Function to submit Vina job for a single control
submit_vina_control() {
    local target=$1
    local control_id=$2
    local smiles=$3
    local target_dir=$4

    local INPUT_CSV="${target_dir}/input/control.csv"
    local OUTPUT_DIR="${target_dir}/Vina/output"
    local TOKEN_FILE="${OUTPUT_DIR}/vina.done"
    local PDB_FILE="${TASK_ROOT}/Input/protein_file/${target}/${target}.pdb"

    # Create input CSV with single SMILES
    echo "ligand_description" > "$INPUT_CSV"
    echo "$smiles" >> "$INPUT_CSV"

    # Check if PDB file exists
    if [ ! -f "$PDB_FILE" ]; then
        log_error "PDB file not found: $PDB_FILE"
        return 1
    fi

    # Get docking box parameters
    BOX_PARAMS=$(get_box_params "$target")
    if [ $? -ne 0 ]; then
        log_error "Failed to get box parameters for ${target}"
        return 1
    fi

    BOX_CENTER=$(echo "$BOX_PARAMS" | head -1)
    BOX_SIZE=$(echo "$BOX_PARAMS" | tail -1)

    if [ -z "$BOX_CENTER" ] || [ -z "$BOX_SIZE" ]; then
        log_error "Failed to extract docking box parameters for ${target}"
        return 1
    fi

    # Create SLURM job script
    local JOB_SCRIPT="${OUTPUT_DIR}/slurm_vina.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=vina_ctrl_${target}_${control_id}
#SBATCH --constraint=${CONSTRAINT}
#SBATCH --time=${VINA_TIME_LIMIT}
#SBATCH --mem=${VINA_MEMORY}
#SBATCH --cpus-per-task=${VINA_CPUS}
#SBATCH --output=${OUTPUT_DIR}/slurm_%j.out
#SBATCH --error=${OUTPUT_DIR}/slurm_%j.err

set -e

echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Control: ${target} control_${control_id}"
echo "Docking box center: ${BOX_CENTER}"
echo "Docking box size: ${BOX_SIZE}"

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate vina

mkdir -p "${OUTPUT_DIR}"

python "${VINA_EXE}" \\
    --smiles "${INPUT_CSV}" \\
    --pdb "${PDB_FILE}" \\
    --box-center ${BOX_CENTER} \\
    --box-size ${BOX_SIZE} \\
    --output "${OUTPUT_DIR}" \\
    --smiles-col "ligand_description"

if [ \$? -eq 0 ]; then
    echo "Vina docking completed successfully"
    touch "${TOKEN_FILE}"
else
    echo "Vina docking failed"
    exit 1
fi

echo "Job completed at: \$(date)"
EOF

    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted Vina for ${target} control_${control_id} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${CONTROLS_OUTPUT_DIR}/submitted_jobs_vina.txt"
        return 0
    else
        log_error "Failed to submit Vina for ${target} control_${control_id}"
        return 1
    fi
}

# Function to submit Boltz2 job for a single control
submit_boltz2_control() {
    local target=$1
    local control_id=$2
    local smiles=$3
    local target_dir=$4

    local INPUT_DIR="${target_dir}/Boltz2/input"
    local OUTPUT_DIR="${target_dir}/Boltz2/output"
    local TOKEN_FILE="${OUTPUT_DIR}/boltz2.done"
    local YAML_FILE="${INPUT_DIR}/control.yaml"

    # Prefolded MSA file from prefold_boltz2
    local MSA_FILE="/home/ubuntu/${target}/boltz2_tmp/boltz_results_${target}/msa/${target}_0.csv"

    mkdir -p "$INPUT_DIR"
    mkdir -p "$OUTPUT_DIR"

    # Check if prefolded MSA exists
    if [ ! -f "$MSA_FILE" ]; then
        log_error "Prefolded MSA not found for ${target}: ${MSA_FILE}"
        log_error "Please run prefold_boltz2 first"
        return 1
    fi

    # Create temporary CSV with single SMILES for gen_boltz_yaml.py
    local TEMP_SMILES_CSV="${INPUT_DIR}/temp_smiles.csv"
    echo "SMILES" > "$TEMP_SMILES_CSV"
    echo "$smiles" >> "$TEMP_SMILES_CSV"

    # Use gen_boltz_yaml.py to create YAML with prefolded MSA
    python3 "${SCRIPT_ROOT}/gen_boltz_yaml.py" \
        --output "$INPUT_DIR" \
        --msa "$MSA_FILE" \
        --smiles-path "$TEMP_SMILES_CSV" \
        --protein-name "$target" \
        --protein-file "$SEQUENCES_CSV"

    if [ $? -ne 0 ]; then
        log_error "Failed to generate Boltz2 YAML for ${target} control_${control_id}"
        return 1
    fi

    # The script creates 0.yaml, rename it to control.yaml
    if [ -f "${INPUT_DIR}/0.yaml" ]; then
        mv "${INPUT_DIR}/0.yaml" "$YAML_FILE"
    else
        log_error "Generated YAML not found: ${INPUT_DIR}/0.yaml"
        return 1
    fi

    # Create SLURM job script
    local JOB_SCRIPT="${OUTPUT_DIR}/slurm_boltz2.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=boltz2_ctrl_${target}_${control_id}
#SBATCH --constraint=${CONSTRAINT}
#SBATCH --time=${BOLTZ2_TIME_LIMIT}
#SBATCH --mem=${BOLTZ2_MEMORY}
#SBATCH --cpus-per-task=${BOLTZ2_CPUS}
#SBATCH --output=${OUTPUT_DIR}/slurm_%j.out
#SBATCH --error=${OUTPUT_DIR}/slurm_%j.err

set -e

echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Control: ${target} control_${control_id}"

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate boltz

if "${BOLTZ2_EXE}" predict "${YAML_FILE}" \\
    --out_dir="${OUTPUT_DIR}" \\
    --override; then
    echo "Boltz2 prediction completed successfully"
    touch "${TOKEN_FILE}"
else
    echo "Boltz2 prediction failed"
    exit 1
fi

echo "Job completed at: \$(date)"
EOF

    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted Boltz2 for ${target} control_${control_id} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${CONTROLS_OUTPUT_DIR}/submitted_jobs_boltz2.txt"
        return 0
    else
        log_error "Failed to submit Boltz2 for ${target} control_${control_id}"
        return 1
    fi
}

# Function to submit AlphaFold3 job for a single control
submit_af3_control() {
    local target=$1
    local control_id=$2
    local smiles=$3
    local target_dir=$4

    local INPUT_DIR="${target_dir}/AF3/input"
    local OUTPUT_DIR="${target_dir}/AF3/output"
    local TOKEN_FILE="${OUTPUT_DIR}/af3.done"
    local JSON_FILE="${INPUT_DIR}/control.json"

    # Prefolded structure with MSA and templates from prefold_af3
    local TARGET_LOWER=$(echo "$target" | tr '[:upper:]' '[:lower:]')
    local PREFOLD_JSON="${TASK_ROOT}/${target}/fine_screening/AF3/prefold/${TARGET_LOWER}/${TARGET_LOWER}_data.json"

    mkdir -p "$INPUT_DIR"
    mkdir -p "$OUTPUT_DIR"

    # Check if prefolded structure exists
    if [ ! -f "$PREFOLD_JSON" ]; then
        log_error "Prefolded AF3 structure not found for ${target}: ${PREFOLD_JSON}"
        log_error "Please run prefold_af3 first"
        return 1
    fi

    # Create temporary CSV with single SMILES for gen_af3_json_with_cmpds.py
    local TEMP_SMILES_CSV="${INPUT_DIR}/temp_smiles.csv"
    echo "SMILES" > "$TEMP_SMILES_CSV"
    echo "$smiles" >> "$TEMP_SMILES_CSV"

    # Use gen_af3_json_with_cmpds.py to create JSON with prefolded data
    python3 "${SCRIPT_ROOT}/gen_af3_json_with_cmpds.py" \
        --output-dir "$INPUT_DIR" \
        --input-json "$PREFOLD_JSON" \
        --smiles-file "$TEMP_SMILES_CSV" \
        --smiles-col "SMILES"

    if [ $? -ne 0 ]; then
        log_error "Failed to generate AF3 JSON for ${target} control_${control_id}"
        return 1
    fi

    # The script creates {target}_0.json, rename it to control.json
    local GENERATED_JSON="${INPUT_DIR}/${target}_0.json"
    if [ -f "$GENERATED_JSON" ]; then
        mv "$GENERATED_JSON" "$JSON_FILE"
    else
        log_error "Generated JSON not found: ${GENERATED_JSON}"
        return 1
    fi

    # Create SLURM job script
    local JOB_SCRIPT="${OUTPUT_DIR}/slurm_af3.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=af3_ctrl_${target}_${control_id}
#SBATCH --constraint=${CONSTRAINT}
#SBATCH --time=${AF3_TIME_LIMIT}
#SBATCH --mem=${AF3_MEMORY}
#SBATCH --cpus-per-task=${AF3_CPUS}
#SBATCH --output=${OUTPUT_DIR}/slurm_%j.out
#SBATCH --error=${OUTPUT_DIR}/slurm_%j.err

set +eu

echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Control: ${target} control_${control_id}"

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate af3

# Run AlphaFold3 prediction
python -u "${AF3_EXE}" \\
    --json_path "${JSON_FILE}" \\
    --model_dir /shared/programs/af3_weights \\
    --db_dir /shared/programs/af3_data \\
    --jackhmmer_binary_path /home/ubuntu/Applications/hmmer/bin/jackhmmer \\
    --hmmalign_binary_path /home/ubuntu/Applications/hmmer/bin/hmmalign \\
    --hmmbuild_binary_path /home/ubuntu/Applications/hmmer/bin/hmmbuild \\
    --hmmsearch_binary_path /home/ubuntu/Applications/hmmer/bin/hmmsearch \\
    --nhmmer_binary_path /home/ubuntu/Applications/hmmer/bin/nhmmer \\
    --output_dir "${OUTPUT_DIR}"

if [ \$? -eq 0 ]; then
    echo "AlphaFold3 completed successfully"
    touch "${TOKEN_FILE}"
    EXIT_CODE=0
else
    echo "AlphaFold3 failed"
    EXIT_CODE=1
fi

echo "Job completed at: \$(date)"
exit \$EXIT_CODE
EOF

    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted AF3 for ${target} control_${control_id} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${CONTROLS_OUTPUT_DIR}/submitted_jobs_af3.txt"
        return 0
    else
        log_error "Failed to submit AF3 for ${target} control_${control_id}"
        return 1
    fi
}

# Function to submit DiffDock job for a single control
submit_diffdock_control() {
    local target=$1
    local control_id=$2
    local smiles=$3
    local target_dir=$4

    local INPUT_CSV="${target_dir}/input/control.csv"
    local OUTPUT_DIR="${target_dir}/DiffDock/output"
    local TOKEN_FILE="${OUTPUT_DIR}/diffdock.done"

    # DiffDock uses the same CSV format as Vina
    if [ ! -f "$INPUT_CSV" ]; then
        echo "ligand_description" > "$INPUT_CSV"
        echo "$smiles" >> "$INPUT_CSV"
    fi

    mkdir -p "$OUTPUT_DIR"

    # Create SLURM job script
    local JOB_SCRIPT="${OUTPUT_DIR}/slurm_diffdock.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=diffdock_ctrl_${target}_${control_id}
#SBATCH --constraint=${CONSTRAINT}
#SBATCH --time=${DIFFDOCK_TIME_LIMIT}
#SBATCH --mem=${DIFFDOCK_MEMORY}
#SBATCH --cpus-per-task=${DIFFDOCK_CPUS}
#SBATCH --output=${OUTPUT_DIR}/slurm_%j.out
#SBATCH --error=${OUTPUT_DIR}/slurm_%j.err

set +eu

echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Control: ${target} control_${control_id}"

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate diffdock

cd "${DIFFDOCK_DIR}"

python -m inference --config "${DIFFDOCK_CONFIG}" \\
    --protein_ligand_csv "${INPUT_CSV}" \\
    --out_dir "${OUTPUT_DIR}"

INFERENCE_EXIT_CODE=\$?

if [ \$INFERENCE_EXIT_CODE -eq 0 ]; then
    echo "DiffDock inference completed successfully"
    touch "${TOKEN_FILE}"
    EXIT_CODE=0
else
    echo "DiffDock inference failed with exit code \$INFERENCE_EXIT_CODE"
    EXIT_CODE=1
fi

echo "Job completed at: \$(date)"
exit \$EXIT_CODE
EOF

    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted DiffDock for ${target} control_${control_id} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${CONTROLS_OUTPUT_DIR}/submitted_jobs_diffdock.txt"
        return 0
    else
        log_error "Failed to submit DiffDock for ${target} control_${control_id}"
        return 1
    fi
}

# Function to submit MD+PBSA job for a single control
submit_md_pbsa_control() {
    local target=$1
    local control_id=$2
    local target_dir=$3

    local DIFFDOCK_OUTPUT="${target_dir}/DiffDock/output"
    local SDF_FILE="${DIFFDOCK_OUTPUT}/${target}_control_${control_id}/rank1.sdf"
    local MD_OUTDIR="${target_dir}/MD_PBSA/MD"
    local MD_TOKEN="${MD_OUTDIR}/token.done"
    local PBSA_OUTDIR="${target_dir}/MD_PBSA/PBSA"
    local PBSA_TOKEN="${PBSA_OUTDIR}/token.done"

    local GRO_FILE="${TASK_ROOT}/Input/protein_file/${target}/${target}.gro"
    local TOP_FILE="${TASK_ROOT}/Input/protein_file/${target}/system_EM.top"

    # Check required files
    if [ ! -f "$GRO_FILE" ]; then
        log_error "GRO file not found: $GRO_FILE"
        return 1
    fi

    if [ ! -f "$TOP_FILE" ]; then
        log_error "TOP file not found: $TOP_FILE"
        return 1
    fi

    mkdir -p "$MD_OUTDIR"
    mkdir -p "$PBSA_OUTDIR"

    # Create SLURM job script
    local JOB_SCRIPT="${target_dir}/MD_PBSA/slurm_md_pbsa.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=md_pbsa_ctrl_${target}_${control_id}
#SBATCH --constraint=g5.16xlarge
#SBATCH --time=${MD_PBSA_TIME_LIMIT}
#SBATCH --cpus-per-task=${MD_PBSA_CPUS}
#SBATCH --output=${target_dir}/MD_PBSA/slurm_%j.out
#SBATCH --error=${target_dir}/MD_PBSA/slurm_%j.err

echo "=========================================="
echo "MD + PBSA Combined Job"
echo "=========================================="
echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Control: ${target} control_${control_id}"
echo ""

# Check if SDF file exists
if [ ! -f "${SDF_FILE}" ]; then
    echo "=========================================="
    echo "SKIPPING JOB: SDF file not found"
    echo "=========================================="
    echo "Expected SDF file: ${SDF_FILE}"
    echo "DiffDock must complete first."
    echo "Skipping MD+PBSA for ${target} control_${control_id}"
    echo "Job completed at: \$(date)"
    echo "=========================================="
    exit 0
fi

set -e

# ============================================================================
# STEP 1: Run Molecular Dynamics
# ============================================================================
echo "=========================================="
echo "STEP 1: Running Molecular Dynamics"
echo "=========================================="

# Create local temporary directory on NVMe/SSD to avoid I/O contention
LOCAL_MD_DIR=/tmp/md_ctrl_${target}_${control_id}_\${SLURM_JOB_ID}
mkdir -p \$LOCAL_MD_DIR

echo "Using local temporary directory: \$LOCAL_MD_DIR"
echo "This avoids disk I/O contention on shared filesystem"
echo ""

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate gmxMMPBSA

bash "${MD_SCRIPT}" \\
    "${SDF_FILE}" \\
    "\$LOCAL_MD_DIR" \\
    "${GRO_FILE}" \\
    "${TOP_FILE}" \\
    "${PBSA_SCRIPT_DIR}"

if [ \$? -eq 0 ]; then
    echo "MD step completed successfully"

    # Copy results from local to network storage
    echo "Copying MD results from local storage to network storage..."
    mkdir -p "${MD_OUTDIR}"
    rsync -avz --progress "\$LOCAL_MD_DIR/" "${MD_OUTDIR}/"

    if [ \$? -eq 0 ]; then
        echo "MD results copied successfully"
        touch "${MD_TOKEN}"
    else
        echo "Failed to copy MD results"
        rm -rf "\$LOCAL_MD_DIR"
        exit 1
    fi
else
    echo "MD step failed"
    rm -rf "\$LOCAL_MD_DIR"
    exit 1
fi

# ============================================================================
# STEP 2: Run PBSA Analysis
# ============================================================================
echo "=========================================="
echo "STEP 2: Running PBSA Analysis"
echo "=========================================="

XTC_FILE="${MD_OUTDIR}/T298.xtc"
TPR_FILE="${MD_OUTDIR}/T298.tpr"
MD_TOP_FILE="${MD_OUTDIR}/system.top"
NDX_FILE="${MD_OUTDIR}/index.ndx"
DAT_FILE="${PBSA_OUTDIR}/FINAL_RESULTS_MMPBSA.dat"

if [ ! -f "\${XTC_FILE}" ] || [ ! -f "\${TPR_FILE}" ] || [ ! -f "\${MD_TOP_FILE}" ] || [ ! -f "\${NDX_FILE}" ]; then
    echo "ERROR: MD output files not found"
    rm -rf "\$LOCAL_MD_DIR"
    exit 1
fi

# Clean up local MD directory
echo "Cleaning up local temporary directory..."
rm -rf "\$LOCAL_MD_DIR"
echo ""

set +eu
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate gmxMMPBSA
bash -c "source ${GMX_RC}"
export PATH="/home/ubuntu/miniconda3/envs/gmxMMPBSA/bin:\$PATH"

cd "${PBSA_OUTDIR}"

env "PATH=\$PATH" "${PBSA_EXE}" -O \\
    -i "${PBSA_SCRIPT_DIR}/mmpbsa.in" \\
    -cs "\${TPR_FILE}" \\
    -ct "\${XTC_FILE}" \\
    -ci "\${NDX_FILE}" \\
    -cg 1 13 \\
    -cp "\${MD_TOP_FILE}" \\
    -o "\${DAT_FILE}" \\
    -eo "${PBSA_OUTDIR}/FINAL_RESULTS_MMPBSA.csv"

if [ \$? -eq 0 ]; then
    echo "PBSA step completed successfully"
    touch "${PBSA_TOKEN}"
    EXIT_CODE=0
else
    echo "PBSA step failed"
    EXIT_CODE=1
fi

echo "=========================================="
echo "Job completed at: \$(date)"
echo "=========================================="
exit \$EXIT_CODE
EOF

    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

    if [ -n "$JOB_ID" ]; then
        log_info "Submitted MD+PBSA for ${target} control_${control_id} (Job ID: ${JOB_ID})"
        echo "$JOB_ID" >> "${CONTROLS_OUTPUT_DIR}/submitted_jobs_md_pbsa.txt"
        return 0
    else
        log_error "Failed to submit MD+PBSA for ${target} control_${control_id}"
        return 1
    fi
}

# ============================================================================
# Prerequisite Checks
# ============================================================================

check_prerequisites() {
    log_info "Checking prerequisites"

    local all_ok=1

    # Check TASK_ROOT exists
    if [ ! -d "$TASK_ROOT" ]; then
        log_error "TASK_ROOT directory not found: $TASK_ROOT"
        all_ok=0
    else
        echo "✓ TASK_ROOT exists: $TASK_ROOT"
    fi

    # Check controls.csv exists
    if [ ! -f "$CONTROLS_CSV" ]; then
        log_error "Controls CSV not found: $CONTROLS_CSV"
        all_ok=0
    else
        echo "✓ Controls CSV found: $CONTROLS_CSV"
        local n_controls=$(tail -n +2 "$CONTROLS_CSV" | wc -l)
        echo "  Found ${n_controls} control molecules"
    fi

    # Check sequences.csv exists
    if [ ! -f "$SEQUENCES_CSV" ]; then
        log_error "Sequences CSV not found: $SEQUENCES_CSV"
        all_ok=0
    else
        echo "✓ Sequences CSV found: $SEQUENCES_CSV"
    fi

    echo ""

    if [ $all_ok -eq 0 ]; then
        log_error "Prerequisites check failed"
        return 1
    fi

    log_info "Prerequisites check passed"
    return 0
}

# ============================================================================
# Help Function
# ============================================================================

show_help() {
    cat <<EOF
Control Molecules Workflow Script
==================================

Usage: $0 [OPTIONS]

This script processes control molecules from task_root/Input/controls.csv
and runs selected workflows (Vina, Boltz2, AF3, DiffDock, MD+PBSA) for each control.

Input File Format (controls.csv):
  Columns: target, SMILES
  - target: Must match a protein name in sequences.csv
  - SMILES: SMILES string of the control molecule

Prerequisites:
  - For Boltz2: Prefolded structures must exist at:
    /home/ubuntu/{target}/boltz2_tmp/boltz_results_{target}/msa/{target}_0.csv
    Run 'prefold_boltz2' from Snakefile first

  - For AF3: Prefolded structures must exist at:
    {task_root}/{target}/fine_screening/AF3/prefold/{target_lower}/{target_lower}_data.json
    Run 'prefold_af3' from Snakefile first

  - For MD+PBSA: DiffDock must complete first to generate SDF files

Options:
  --task-root DIR      Set task root directory (default: $TASK_ROOT)
  --skip-vina          Skip Vina workflow
  --skip-boltz2        Skip Boltz2 workflow
  --skip-af3           Skip AlphaFold3 workflow
  --skip-diffdock      Skip DiffDock workflow
  --skip-md-pbsa       Skip MD+PBSA workflow
  --help               Show this help message

Examples:
  # Run all workflows (default)
  $0

  # Run only Vina and Boltz2
  $0 --skip-af3 --skip-diffdock --skip-md-pbsa

  # Use custom task root
  $0 --task-root /path/to/data

Configuration:
  Task Root: $TASK_ROOT
  Controls CSV: $CONTROLS_CSV
  Sequences CSV: $SEQUENCES_CSV

  Workflows enabled:
    Vina: $([ $RUN_VINA -eq 1 ] && echo "Yes" || echo "No")
    Boltz2: $([ $RUN_BOLTZ2 -eq 1 ] && echo "Yes" || echo "No")
    AlphaFold3: $([ $RUN_AF3 -eq 1 ] && echo "Yes" || echo "No")
    DiffDock: $([ $RUN_DIFFDOCK -eq 1 ] && echo "Yes" || echo "No")
    MD+PBSA: $([ $RUN_MD_PBSA -eq 1 ] && echo "Yes" || echo "No")

Output Structure:
  ${TASK_ROOT}/Controls/
    ├── <target>/
    │   └── control_<id>/
    │       ├── input/
    │       ├── Vina/output/
    │       ├── Boltz2/{input,output}/
    │       ├── AF3/{input,output}/
    │       ├── DiffDock/output/
    │       └── MD_PBSA/{MD,PBSA}/

Input File Generation:
  - Boltz2: Uses scripts/gen_boltz_yaml.py with prefolded MSA
  - AF3: Uses scripts/gen_af3_json_with_cmpds.py with prefolded data.json

Notes:
  - MD+PBSA requires DiffDock to complete first (produces SDF files)
  - Boltz2 and AF3 use prefolded protein structures with MSA and templates
  - All algorithms run by default, use --skip-* flags to disable
  - Job IDs are saved to ${CONTROLS_OUTPUT_DIR}/submitted_jobs_*.txt

EOF
}

# ============================================================================
# Main Script
# ============================================================================

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --task-root)
            TASK_ROOT="$2"
            CONTROLS_CSV="${TASK_ROOT}/Input/controls.csv"
            SEQUENCES_CSV="${TASK_ROOT}/Input/sequences.csv"
            CONTROLS_OUTPUT_DIR="${TASK_ROOT}/Controls"
            shift 2
            ;;
        --skip-vina)
            RUN_VINA=0
            shift
            ;;
        --skip-boltz2)
            RUN_BOLTZ2=0
            shift
            ;;
        --skip-af3)
            RUN_AF3=0
            shift
            ;;
        --skip-diffdock)
            RUN_DIFFDOCK=0
            shift
            ;;
        --skip-md-pbsa)
            RUN_MD_PBSA=0
            shift
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run '$0 --help' for usage information"
            exit 1
            ;;
    esac
done

# Show banner
echo "=========================================="
echo "  Control Molecules Workflow"
echo "=========================================="
echo ""
echo "Task Root: $TASK_ROOT"
echo "Controls CSV: $CONTROLS_CSV"
echo ""

# Check prerequisites
if ! check_prerequisites; then
    log_error "Prerequisites check failed. Please fix the issues and try again."
    exit 1
fi

# Create main output directory
mkdir -p "$CONTROLS_OUTPUT_DIR"

# Clear previous job lists
> "${CONTROLS_OUTPUT_DIR}/submitted_jobs_vina.txt"
> "${CONTROLS_OUTPUT_DIR}/submitted_jobs_boltz2.txt"
> "${CONTROLS_OUTPUT_DIR}/submitted_jobs_af3.txt"
> "${CONTROLS_OUTPUT_DIR}/submitted_jobs_diffdock.txt"
> "${CONTROLS_OUTPUT_DIR}/submitted_jobs_md_pbsa.txt"

# Read controls CSV and process each control
log_info "Processing control molecules"

CONTROL_ID=0
while IFS=, read -r target smiles; do
    # Skip header
    if [ "$target" == "target" ]; then
        continue
    fi

    log_info "Processing control ${CONTROL_ID}: target=${target}, SMILES=${smiles}"

    # Setup directories
    TARGET_DIR=$(setup_control_dirs "$target" "$CONTROL_ID")

    # Submit workflows
    if [ $RUN_VINA -eq 1 ]; then
        submit_vina_control "$target" "$CONTROL_ID" "$smiles" "$TARGET_DIR" || true
    fi

    if [ $RUN_BOLTZ2 -eq 1 ]; then
        submit_boltz2_control "$target" "$CONTROL_ID" "$smiles" "$TARGET_DIR" || true
    fi

    if [ $RUN_AF3 -eq 1 ]; then
        submit_af3_control "$target" "$CONTROL_ID" "$smiles" "$TARGET_DIR" || true
    fi

    if [ $RUN_DIFFDOCK -eq 1 ]; then
        submit_diffdock_control "$target" "$CONTROL_ID" "$smiles" "$TARGET_DIR" || true
    fi

    if [ $RUN_MD_PBSA -eq 1 ]; then
        submit_md_pbsa_control "$target" "$CONTROL_ID" "$TARGET_DIR" || true
    fi

    CONTROL_ID=$((CONTROL_ID + 1))

done < "$CONTROLS_CSV"

# Summary
log_info "Control workflow submission completed!"
echo "Processed $CONTROL_ID control molecules"
echo ""
echo "Job IDs saved to:"
echo "  Vina:     ${CONTROLS_OUTPUT_DIR}/submitted_jobs_vina.txt"
echo "  Boltz2:   ${CONTROLS_OUTPUT_DIR}/submitted_jobs_boltz2.txt"
echo "  AF3:      ${CONTROLS_OUTPUT_DIR}/submitted_jobs_af3.txt"
echo "  DiffDock: ${CONTROLS_OUTPUT_DIR}/submitted_jobs_diffdock.txt"
echo "  MD+PBSA:  ${CONTROLS_OUTPUT_DIR}/submitted_jobs_md_pbsa.txt"
echo ""
echo "Monitor jobs with: squeue -u \$USER"
echo "Check outputs in: ${CONTROLS_OUTPUT_DIR}/"
echo ""
