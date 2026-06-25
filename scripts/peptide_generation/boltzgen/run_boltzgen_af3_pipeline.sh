#!/bin/bash
# BoltzGen → AF3 peptide design pipeline
# Step 1: BoltzGen design (replaces RFDiffusion + ProteinMPNN)
# Steps 3-6: shared with the RFDiffusion pipeline

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Defaults
PROTEIN=""
PROTEIN_PDB=""
PEPTIDE_LENGTH=""
OUTPUT_DIR=""
HOTSPOT_RESIDUES=""
N_DESIGNS=200
N_BATCHES=40
RESUME_FROM=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --protein)           PROTEIN="$2";           shift 2 ;;
        --protein-pdb)       PROTEIN_PDB="$2";       shift 2 ;;
        --peptide-length)    PEPTIDE_LENGTH="$2";    shift 2 ;;
        --output-dir)        OUTPUT_DIR="$2";        shift 2 ;;
        --hotspot-residues)  HOTSPOT_RESIDUES="$2";  shift 2 ;;
        --n-designs)         N_DESIGNS="$2";         shift 2 ;;
        --n-batches)         N_BATCHES="$2";         shift 2 ;;
        --resume-from)       RESUME_FROM="$2";       shift 2 ;;
        -h|--help)
            cat <<HELP
Usage: $0 [OPTIONS]

Required:
  --protein NAME              Protein name
  --protein-pdb PATH          Path to protein PDB or CIF structure
  --peptide-length N          Target peptide length (passed to BoltzGen as N..N)
  --output-dir PATH           Output directory

Optional:
  --hotspot-residues RES      Hotspot residues e.g. "A:45,A:67,A:123"
  --n-designs N               Number of BoltzGen designs (default: 200)
  --n-batches N               AF3 SLURM batch count (default: 40)
  --resume-from STEP          Resume from step number (1-6)

Example:
  $0 --protein JAK1 \\
     --protein-pdb /path/to/jak1.pdb \\
     --peptide-length 15 \\
     --output-dir /path/to/output \\
     --hotspot-residues "A:45,A:67"
HELP
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$PROTEIN" || -z "$PROTEIN_PDB" || -z "$PEPTIDE_LENGTH" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required arguments. Use --help for usage."
    exit 1
fi

if [[ ! -f "$PROTEIN_PDB" ]]; then
    echo "Error: Protein structure not found: $PROTEIN_PDB"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
STATE_FILE="${OUTPUT_DIR}/.pipeline_state.json"

if [[ ! -f "$STATE_FILE" ]]; then
    cat > "$STATE_FILE" <<JSON
{
    "protein": "$PROTEIN",
    "protein_pdb": "$PROTEIN_PDB",
    "peptide_length": $PEPTIDE_LENGTH,
    "n_designs": $N_DESIGNS,
    "hotspot_residues": "$HOTSPOT_RESIDUES",
    "completed_steps": []
}
JSON
fi

log_info() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log_step() {
    echo ""
    echo "=========================================="
    echo "STEP $1: $2"
    echo "=========================================="
}

mark_complete() {
    jq --arg step "$1" '.completed_steps += [$step] | .completed_steps |= unique' \
        "$STATE_FILE" > "${STATE_FILE}.tmp"
    mv "${STATE_FILE}.tmp" "$STATE_FILE"
}

is_complete() {
    jq -e --arg step "$1" '.completed_steps | index($step)' "$STATE_FILE" > /dev/null 2>&1
}

extract_protein_fasta() {
    local pdb_file=$1
    local output_fasta=$2

    python3 <<PYTHON_EXTRACT
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder
import os

path = '$pdb_file'
ext = os.path.splitext(path)[1].lower()
if ext in ('.cif', '.mmcif'):
    parser = MMCIFParser(QUIET=True)
else:
    parser = PDBParser(QUIET=True)

structure = parser.get_structure('protein', path)
ppb = PPBuilder()

records = []
for chain in structure[0]:
    seq = ''.join(str(pp.get_sequence()) for pp in ppb.build_peptides(chain))
    if seq:
        records.append((chain.id, seq))

with open('$output_fasta', 'w') as f:
    for chain_id, seq in records:
        f.write(f'>${PROTEIN}_{chain_id}\n')
        f.write(seq + '\n')

print(f"Extracted {len(records)} chain(s) to $output_fasta")
PYTHON_EXTRACT
}

log_info "Starting BoltzGen → AF3 pipeline for $PROTEIN"
log_info "Output directory: $OUTPUT_DIR"
log_info "Parameters: ${N_DESIGNS} designs, peptide length ${PEPTIDE_LENGTH}"

# Step 1: BoltzGen design
if [[ $RESUME_FROM -le 1 ]] && ! is_complete "step1"; then
    log_step 1 "BoltzGen Peptide Design"

    bash "${SCRIPT_DIR}/step1_run_boltzgen.sh" \
        --protein-pdb    "$PROTEIN_PDB" \
        --output-dir     "$OUTPUT_DIR" \
        --peptide-length "$PEPTIDE_LENGTH" \
        --n-designs      "$N_DESIGNS" \
        ${HOTSPOT_RESIDUES:+--hotspot-residues "$HOTSPOT_RESIDUES"}

    mark_complete "step1"
else
    log_info "Skipping Step 1 (already complete or resume point)"
fi

# Step 2: Extract sequences from BoltzGen output
if [[ $RESUME_FROM -le 2 ]] && ! is_complete "step2"; then
    log_step 2 "Extract Sequences from BoltzGen Output"

    python3 "${SCRIPT_DIR}/step2_extract_sequences.py" \
        --boltzgen-output "${OUTPUT_DIR}/01_boltzgen/output" \
        --output-dir      "${OUTPUT_DIR}/02_sequences"

    mark_complete "step2"
else
    log_info "Skipping Step 2 (already complete or resume point)"
fi

# Step 3: Prepare AF3 inputs
if [[ $RESUME_FROM -le 3 ]] && ! is_complete "step3"; then
    log_step 3 "Prepare AF3 Input Files"

    PROTEIN_FASTA="${OUTPUT_DIR}/protein.fasta"
    extract_protein_fasta "$PROTEIN_PDB" "$PROTEIN_FASTA"

    python3 "${PARENT_DIR}/step3_prepare_af3.py" \
        --protein-fasta  "$PROTEIN_FASTA" \
        --peptide-dir    "${OUTPUT_DIR}/02_sequences/seqs" \
        --output-dir     "${OUTPUT_DIR}/03_af3_input" \
        --protein-name   "$PROTEIN" \
        --no-skip-native

    mark_complete "step3"
else
    log_info "Skipping Step 3 (already complete or resume point)"
fi

# Step 4: Submit AF3 SLURM jobs
if [[ $RESUME_FROM -le 4 ]] && ! is_complete "step4"; then
    log_step 4 "Submit AF3 Validation Jobs"

    bash "${PARENT_DIR}/step4_submit_af3.sh" \
        --output-dir    "$OUTPUT_DIR" \
        --protein-name  "$PROTEIN" \
        --n-batches     "$N_BATCHES"

    mark_complete "step4"
else
    log_info "Skipping Step 4 (already complete or resume point)"
fi

# Step 5: Stream AF3 results
if [[ $RESUME_FROM -le 5 ]] && ! is_complete "step5"; then
    log_step 5 "Stream AF3 Results Collection"

    bash "${PARENT_DIR}/step5_collect_streaming.sh" \
        --output-dir    "$OUTPUT_DIR" \
        --protein-name  "$PROTEIN" \
        --poll-interval 60 \
        --max-iterations 1000

    mark_complete "step5"
else
    log_info "Skipping Step 5 (already complete or resume point)"
fi

# Step 6: Finalize ranking
if [[ $RESUME_FROM -le 6 ]] && ! is_complete "step6"; then
    log_step 6 "Finalize Peptide Ranking"

    python3 "${PARENT_DIR}/step6_finalize_ranking.py" \
        --output-dir    "$OUTPUT_DIR" \
        --protein-name  "$PROTEIN"

    mark_complete "step6"
else
    log_info "Skipping Step 6 (already complete or resume point)"
fi

log_info ""
log_info "=========================================="
log_info "Pipeline Complete!"
log_info "=========================================="
log_info "Results: ${OUTPUT_DIR}/05_results/ranked_peptides.csv"
