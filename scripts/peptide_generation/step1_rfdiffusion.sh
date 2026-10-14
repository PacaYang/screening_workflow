#!/bin/bash
# Step 1: RFDiffusion backbone generation

set -e

# Parse arguments
PROTEIN_PDB=""
OUTPUT_DIR=""
PEPTIDE_LENGTH=""
N_DESIGNS=50
HOTSPOT_RESIDUES=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --protein-pdb) PROTEIN_PDB="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --peptide-length) PEPTIDE_LENGTH="$2"; shift 2 ;;
        --n-designs) N_DESIGNS="$2"; shift 2 ;;
        --hotspot-residues) HOTSPOT_RESIDUES="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Validate inputs
if [[ -z "$PROTEIN_PDB" || -z "$OUTPUT_DIR" || -z "$PEPTIDE_LENGTH" ]]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --protein-pdb PDB --output-dir DIR --peptide-length N [--n-designs N] [--hotspot-residues RES]"
    exit 1
fi

if [[ ! -f "$PROTEIN_PDB" ]]; then
    echo "Error: Protein PDB not found: $PROTEIN_PDB"
    exit 1
fi

# Setup paths
RFD_DIR="${OUTPUT_DIR}/01_rfdiffusion"
BACKBONE_DIR="${RFD_DIR}/backbones"
LOG_DIR="${RFD_DIR}/logs"

mkdir -p "$BACKBONE_DIR" "$LOG_DIR"

# Build contig map covering ALL target chains, then the designed peptide.
# Each contiguous run of resolved residues in a chain becomes one
# "<chain><start>-<end>/0" segment (the /0 is a chain break), and the final
# "<len>-<len>" segment is the de novo peptide (gets the next chain letter,
# e.g. D for an A/B/C trimer target).
#
# We split each chain at gaps (missing/unresolved residues) because RFDiffusion
# requires every residue in a contig range to be present in the PDB; a single
# "<first>-<last>" span across a gap makes it assert on the missing residue.
TARGET_SEGMENTS=""
ALL_CHAINS=$(grep "^ATOM" "$PROTEIN_PDB" | awk '{print $5}' | awk '!seen[$0]++')
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Target chains detected: $(echo $ALL_CHAINS | tr '\n' ' ')"

for ch in $ALL_CHAINS; do
    # Emit one "<chain><start>-<end>/0" token per contiguous run of residues.
    ch_runs=$(grep "^ATOM" "$PROTEIN_PDB" | awk -v c="$ch" '$5==c {print $6}' \
        | awk '!seen[$0]++' | sort -n \
        | awk -v ch="$ch" '
            NR==1              {start=$1; prev=$1; next}
            $1==prev+1         {prev=$1; next}
                               {printf "%s%d-%d/0 ", ch, start, prev; start=$1; prev=$1}
            END                {printf "%s%d-%d/0 ", ch, start, prev}')
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]   chain ${ch}: ${ch_runs}"
    TARGET_SEGMENTS="${TARGET_SEGMENTS}${ch_runs}"
done

# Build contig map: all target chains followed by the de novo peptide segment
CONTIG="${TARGET_SEGMENTS}${PEPTIDE_LENGTH}-${PEPTIDE_LENGTH}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Contig map: ${CONTIG}"

# Activate RFDiffusion environment
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Activating conda environment SE3nv..."
set +u  # Temporarily disable unbound variable check for conda activation
source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate SE3nv
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Conda environment activated: $(which python)"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting RFDiffusion with ${N_DESIGNS} designs..."

# Run RFDiffusion
cd ~/Applications/RFdiffusion

if [[ -n "$HOTSPOT_RESIDUES" ]]; then
    # Convert "A:45,A:67" to "A45,A67" format
    HOTSPOT_RFD=$(echo "$HOTSPOT_RESIDUES" | sed 's/://g')
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using hotspot residues: ${HOTSPOT_RFD}"

    python scripts/run_inference.py \
        inference.output_prefix="${BACKBONE_DIR}/design" \
        inference.input_pdb="$PROTEIN_PDB" \
        "contigmap.contigs=[${CONTIG}]" \
        inference.num_designs="$N_DESIGNS" \
        "ppi.hotspot_res=[${HOTSPOT_RFD}]" \
        2>&1 | tee "${LOG_DIR}/rfdiffusion.log"
else
    python scripts/run_inference.py \
        inference.output_prefix="${BACKBONE_DIR}/design" \
        inference.input_pdb="$PROTEIN_PDB" \
        "contigmap.contigs=[${CONTIG}]" \
        inference.num_designs="$N_DESIGNS" \
        2>&1 | tee "${LOG_DIR}/rfdiffusion.log"
fi

# Count generated backbones
N_GENERATED=$(find "$BACKBONE_DIR" -name "design_*.pdb" | wc -l)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generated ${N_GENERATED} backbones in ${BACKBONE_DIR}"

if [[ $N_GENERATED -eq 0 ]]; then
    echo "Error: No backbones generated"
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] RFDiffusion completed successfully"
