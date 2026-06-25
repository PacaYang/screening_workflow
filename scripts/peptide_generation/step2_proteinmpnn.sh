#!/bin/bash
# Step 2: ProteinMPNN sequence design

set -euo pipefail

# Parse arguments
OUTPUT_DIR=""
N_SEQS=10

while [[ $# -gt 0 ]]; do
    case $1 in
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --n-seqs) N_SEQS="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing --output-dir"
    exit 1
fi

# Setup paths
RFD_DIR="${OUTPUT_DIR}/01_rfdiffusion"
BACKBONE_DIR="${RFD_DIR}/backbones"
MPNN_DIR="${OUTPUT_DIR}/02_proteinmpnn"
SEQ_DIR="${MPNN_DIR}/sequences"
LOG_DIR="${MPNN_DIR}/logs"

mkdir -p "$SEQ_DIR" "$LOG_DIR"

if [[ ! -d "$BACKBONE_DIR" ]]; then
    echo "Error: Backbone directory not found: $BACKBONE_DIR"
    exit 1
fi

# Activate ProteinMPNN environment (use SE3nv which has ProteinMPNN)
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
set +u
conda activate SE3nv
set -u

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Parsing PDB chains..."

# Parse chains
cd ~/Applications/ProteinMPNN
python helper_scripts/parse_multiple_chains.py \
    --input_path="$BACKBONE_DIR" \
    --output_path="${MPNN_DIR}/parsed_pdbs.jsonl" \
    2>&1 | tee "${LOG_DIR}/parse_chains.log"

# Determine chain layout from the first backbone. RFDiffusion writes the
# de novo peptide as its own chain (e.g. D for an A/B/C trimer target), but the
# chain ORDER in the PDB is not reliable (RFDiffusion often lists the peptide
# first). The peptide is always the SHORTEST chain, so we identify the design
# chain by smallest residue (CA) count and fix all the longer target chains.
FIRST_PDB=$(find "$BACKBONE_DIR" -name "design_*.pdb" | sort | head -1)

# Build "chain count" pairs (CA atoms per chain), then pick min as peptide.
DESIGN_CHAIN=$(grep "^ATOM" "$FIRST_PDB" | awk '$3=="CA"{c[$5]++} END{for (ch in c) print c[ch], ch}' | sort -n | head -1 | awk '{print $2}')
ALL_CHAINS=$(grep "^ATOM" "$FIRST_PDB" | awk '{print $5}' | awk '!seen[$0]++' | sort)
FIXED_CHAINS=$(echo "$ALL_CHAINS" | grep -v "^${DESIGN_CHAIN}$" | tr '\n' ' ' | sed 's/ *$//')
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Chains: $(echo $ALL_CHAINS | tr '\n' ' ')"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Per-chain CA counts:"
grep "^ATOM" "$FIRST_PDB" | awk '$3=="CA"{c[$5]++} END{for (ch in c) printf "  chain %s: %d CA\n", ch, c[ch]}' | sort
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fixed (target): ${FIXED_CHAINS} | Design (peptide): ${DESIGN_CHAIN}"

# Assign which chain ProteinMPNN is allowed to redesign (peptide only).
# This produces the chain_id_jsonl: target chains stay fixed at their input
# sequence, only the peptide chain is sampled.
python helper_scripts/assign_fixed_chains.py \
    --input_path="${MPNN_DIR}/parsed_pdbs.jsonl" \
    --output_path="${MPNN_DIR}/assigned_chains.jsonl" \
    --chain_list="$DESIGN_CHAIN" \
    2>&1 | tee "${LOG_DIR}/assign_chains.log"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running ProteinMPNN (${N_SEQS} seqs per backbone)..."

# Run ProteinMPNN: only the peptide chain is designed, target chains fixed
python protein_mpnn_run.py \
    --jsonl_path="${MPNN_DIR}/parsed_pdbs.jsonl" \
    --chain_id_jsonl="${MPNN_DIR}/assigned_chains.jsonl" \
    --out_folder="$SEQ_DIR" \
    --num_seq_per_target="$N_SEQS" \
    --sampling_temp="0.1" \
    --seed=37 \
    --batch_size=1 \
    --save_score=1 \
    --save_probs=0 \
    --omit_AAs='X' \
    2>&1 | tee "${LOG_DIR}/proteinmpnn.log"

# Count generated sequences
N_FASTA=$(find "$SEQ_DIR" -name "*.fa" | wc -l)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generated ${N_FASTA} FASTA files in ${SEQ_DIR}"

if [[ $N_FASTA -eq 0 ]]; then
    echo "Error: No sequences generated"
    exit 1
fi

# Aggregate scores
if [ -f "${SEQ_DIR}/seqs/design_0001.fa" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Aggregating ProteinMPNN scores..."

    echo "backbone_id,seq_idx,score,seq_recovery" > "${MPNN_DIR}/scores.csv"

    for score_file in "${SEQ_DIR}"/seqs/*.fa; do
        if grep -q "score" "$score_file"; then
            basename=$(basename "$score_file" .fa)
            grep "score" "$score_file" | awk -v bn="$basename" -F'[=,]' '{print bn","NR","$2","$4}' >> "${MPNN_DIR}/scores.csv"
        fi
    done
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ProteinMPNN completed successfully"
