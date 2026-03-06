#!/bin/bash
#SBATCH --job-name=rfaa_test_single
#SBATCH --time=48:00:00
#SBATCH --mem=15G
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --partition=g24
#SBATCH --output=/home/yangl_pacagen_com/Projects/B3/IL6/IL6RA/fine_screening/RoseTTAFold/protein_ligand/test_single/logs/slurm_%j.out
#SBATCH --error=/home/yangl_pacagen_com/Projects/B3/IL6/IL6RA/fine_screening/RoseTTAFold/protein_ligand/test_single/logs/slurm_%j.err

set -e

echo "Job started at: $(date)"
echo "Testing RoseTTAFold single compound prediction for IL6RA"

# ============================================================================
# Configuration
# ============================================================================

PROTEIN="IL6RA"
RFAA_ROOT="/home/yangl_pacagen_com/Applications/RoseTTAFold-All-Atom"
PROJECT_DIR="/home/yangl_pacagen_com/Projects/B3/IL6/IL6RA"
FOLD_OUTPUT="${PROJECT_DIR}/fine_screening/RoseTTAFold/protein_folding/output"
FASTA_FILE="${PROJECT_DIR}/fine_screening/RoseTTAFold/protein_folding/input/IL6RA.fasta"
TEST_DIR="${PROJECT_DIR}/fine_screening/RoseTTAFold/protein_ligand/test_single"
CONFIG_DIR="${TEST_DIR}/config"
LOG_DIR="${TEST_DIR}/logs"

# Test compound (first row of selected.csv)
SMILES="O=C1CC(O)C(O)c2cc(O)cc(O)c21"
COMPOUND_ID="compound_0"

# Database paths
export DB_UR30="/home/yangl_pacagen_com/Applications/model_weights/rosetta_db/UniRef30_2020_06/UniRef30_2020_06"
export DB_BFD="/home/yangl_pacagen_com/Applications/model_weights/rosetta_db/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"

# ============================================================================
# Setup
# ============================================================================

# Activate environment
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate RFAA_test

# Create directories
mkdir -p "$CONFIG_DIR" "$LOG_DIR" "${TEST_DIR}/${COMPOUND_ID}"

# ============================================================================
# Check protein folding prerequisites
# ============================================================================

PROTEIN_MSA_DIR="${FOLD_OUTPUT}/${PROTEIN}_fold/A"

# Only require MSA file - templates are optional due to PSIPRED buffer overflow issue
if [ ! -f "${PROTEIN_MSA_DIR}/t000_.msa0.a3m" ]; then
    echo "ERROR: Missing MSA file: ${PROTEIN_MSA_DIR}/t000_.msa0.a3m"
    echo "Please ensure protein folding has generated the MSA."
    exit 1
fi

echo "MSA file found. Template files (ss2, hhr, atab) are optional and will be skipped if missing."

# ============================================================================
# Copy pre-computed MSA files
# ============================================================================

COMPOUND_MSA_DIR="${TEST_DIR}/${COMPOUND_ID}/${PROTEIN}_ligand_${COMPOUND_ID}/A"
mkdir -p "$COMPOUND_MSA_DIR"

echo "Copying pre-computed MSA files from protein folding output..."
cp "${PROTEIN_MSA_DIR}/t000_.msa0.a3m" "$COMPOUND_MSA_DIR/"
# Copy template files only if they exist
[ -f "${PROTEIN_MSA_DIR}/t000_.ss2" ] && cp "${PROTEIN_MSA_DIR}/t000_.ss2" "$COMPOUND_MSA_DIR/" || echo "  Skipping t000_.ss2 (not found)"
[ -f "${PROTEIN_MSA_DIR}/t000_.hhr" ] && cp "${PROTEIN_MSA_DIR}/t000_.hhr" "$COMPOUND_MSA_DIR/" || echo "  Skipping t000_.hhr (not found)"
[ -f "${PROTEIN_MSA_DIR}/t000_.atab" ] && cp "${PROTEIN_MSA_DIR}/t000_.atab" "$COMPOUND_MSA_DIR/" || echo "  Skipping t000_.atab (not found)"
echo "MSA files copied to: $COMPOUND_MSA_DIR"

# ============================================================================
# Write YAML config
# ============================================================================

CONFIG_FILE="${CONFIG_DIR}/config_${COMPOUND_ID}.yaml"

cat > "$CONFIG_FILE" <<EOF
defaults:
  - base

job_name: "${PROTEIN}_ligand_${COMPOUND_ID}"
output_path: "${TEST_DIR}/${COMPOUND_ID}"

database_params:
  hhdb: "/home/yangl_pacagen_com/Applications/model_weights/rosetta_db/pdb100_2021Mar03/pdb100_2021Mar03"

protein_inputs:
  A:
    fasta_file: "${FASTA_FILE}"

sm_inputs:
  B:
    input: "${SMILES}"
    input_type: "smiles"
EOF

echo "Config written to: $CONFIG_FILE"
echo "Contents:"
cat "$CONFIG_FILE"

# ============================================================================
# Run RFAA inference
# ============================================================================

cd "$RFAA_ROOT"

echo ""
echo "Starting RFAA inference..."
echo "Config dir: ${CONFIG_DIR}"
echo "Config name: config_${COMPOUND_ID}"

python -m rf2aa.run_inference \
    --config-dir "${CONFIG_DIR}" \
    --config-name "config_${COMPOUND_ID}" \
    2>&1 | tee "${LOG_DIR}/inference_${COMPOUND_ID}.log"

echo ""
echo "Inference completed at: $(date)"

# ============================================================================
# Verify output
# ============================================================================

AUX_FILE=$(find "${TEST_DIR}/${COMPOUND_ID}" -name "*_aux.pt" 2>/dev/null | head -1)
if [ -n "$AUX_FILE" ]; then
    echo "SUCCESS: Output file found: $AUX_FILE"
else
    echo "WARNING: No *_aux.pt file found in ${TEST_DIR}/${COMPOUND_ID}"
    echo "Check logs for errors."
    exit 1
fi

echo "Job completed successfully at: $(date)"
