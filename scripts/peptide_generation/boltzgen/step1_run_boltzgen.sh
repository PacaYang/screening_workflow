#!/bin/bash
# Step 1: Generate design spec YAML and submit BoltzGen SLURM job

set -euo pipefail

PROTEIN_PDB=""
OUTPUT_DIR=""
PEPTIDE_LENGTH=""
N_DESIGNS=200
HOTSPOT_RESIDUES=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --protein-pdb)      PROTEIN_PDB="$2";      shift 2 ;;
        --output-dir)       OUTPUT_DIR="$2";       shift 2 ;;
        --peptide-length)   PEPTIDE_LENGTH="$2";   shift 2 ;;
        --n-designs)        N_DESIGNS="$2";        shift 2 ;;
        --hotspot-residues) HOTSPOT_RESIDUES="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$PROTEIN_PDB" || -z "$OUTPUT_DIR" || -z "$PEPTIDE_LENGTH" ]]; then
    echo "Error: Missing required arguments"
    exit 1
fi

BOLTZGEN_DIR="${OUTPUT_DIR}/01_boltzgen"
mkdir -p "${BOLTZGEN_DIR}/output"

SPEC_FILE="${BOLTZGEN_DIR}/design_spec.yaml"
BUDGET=$(( N_DESIGNS < 50 ? N_DESIGNS : 50 ))

# Convert hotspot residues from "A:45,A:67,B:12" to per-chain binding maps.
#
# The residue numbers we receive are PDB *author* numbers (e.g. C:256), but
# BoltzGen indexes residues 1..N by their position in the chain's SEQRES /
# full_sequence (matching gemmi's align_sequence_to_polymer, which is exactly
# what boltzgen.data.parse.pdb_parser.parse_pdb does). For structures whose
# author numbering is offset or has gaps, the raw author number overruns the
# chain length and BoltzGen raises "Specified end ... is higher than the length
# of the chain." So we translate each author number to its 1-indexed sequence
# position here, using the same alignment BoltzGen uses.
build_binding_yaml() {
    local hotspots="$1"
    local pdb_path="$2"
    /home/yangl_pacagen_com/miniconda3/envs/boltzgen/bin/python <<PYTHON
hotspots = "$hotspots"
pdb_path = "$pdb_path"
from collections import defaultdict, OrderedDict
import sys
import gemmi

# Parse "C:256,C:257,B:12" -> {chain: [author_num, ...]} preserving order.
requested = OrderedDict()
for token in hotspots.split(','):
    token = token.strip()
    if not token:
        continue
    if ':' in token:
        chain_id, resnum = token.split(':', 1)
        chain_id = chain_id.strip()
    else:
        chain_id, resnum = 'A', token
    requested.setdefault(chain_id, []).append(int(resnum.strip()))

# Build author-number -> 1-indexed sequence-position map per chain, mirroring
# boltzgen.data.parse.pdb_parser.parse_pdb.
st = gemmi.read_structure(pdb_path)
st.setup_entities()
author_to_seqpos = defaultdict(dict)  # chain_name -> {author_num: label_seq}
for entity in st.entities:
    if entity.entity_type.name != "Polymer":
        continue
    sc_name = entity.subchains[0]
    sc = st[0].get_subchain(sc_name)
    if len(entity.full_sequence) < 1:
        entity.full_sequence = [res.name for res in sc]
    align = gemmi.align_sequence_to_polymer(
        entity.full_sequence, sc, entity.polymer_type, gemmi.AlignmentScoring()
    ).match_string
    chain_name = next((ch.name for ch in st[0] if ch[0].subchain == sc_name), None)
    i = 0
    for j, a in enumerate(align):
        if a == "|":
            author_to_seqpos[chain_name][sc[i].seqid.num] = j + 1
            i += 1

print("      binding_types:")
for chain_id, author_nums in requested.items():
    seqpos = []
    for a in author_nums:
        pos = author_to_seqpos.get(chain_id, {}).get(a)
        if pos is None:
            sys.exit(
                f"Hotspot residue {chain_id}:{a} not found among resolved "
                f"residues of chain {chain_id} in {pdb_path}"
            )
        seqpos.append(pos)
    print(f"        - chain:")
    print(f"            id: {chain_id}")
    print(f"            binding: {','.join(map(str, seqpos))}")
PYTHON
}

# Write the design spec YAML
{
    echo "entities:"
    echo "  - protein:"
    echo "      id: G"
    echo "      sequence: ${PEPTIDE_LENGTH}..${PEPTIDE_LENGTH}"
    echo ""
    echo "  - file:"
    echo "      path: $(realpath "$PROTEIN_PDB")"
    echo "      include: all"
    if [[ -n "$HOTSPOT_RESIDUES" ]]; then
        build_binding_yaml "$HOTSPOT_RESIDUES" "$(realpath "$PROTEIN_PDB")"
    fi
    echo "      structure_groups: \"all\""
} > "$SPEC_FILE"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Design spec written to $SPEC_FILE"
cat "$SPEC_FILE"

LOG_OUT="${BOLTZGEN_DIR}/slurm_%j.out"
LOG_ERR="${BOLTZGEN_DIR}/slurm_%j.err"
DONE_TOKEN="${BOLTZGEN_DIR}/boltzgen.done"
SLURM_SCRIPT="${BOLTZGEN_DIR}/run_boltzgen.slurm"

BOLTZGEN_BIN="/home/yangl_pacagen_com/miniconda3/envs/boltzgen/bin/boltzgen"

# Write SLURM script to file (avoids heredoc escaping issues)
cat > "$SLURM_SCRIPT" <<SLURM_SCRIPT_EOF
#!/bin/bash
#SBATCH --job-name=boltzgen_design
#SBATCH --output=${LOG_OUT}
#SBATCH --error=${LOG_ERR}
#SBATCH --time=24:00:00
#SBATCH --mem=40G
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --partition=g212

set -e
echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"

${BOLTZGEN_BIN} run ${SPEC_FILE} \
    --output ${BOLTZGEN_DIR}/output \
    --protocol peptide-anything \
    --num_designs ${N_DESIGNS} \
    --budget ${BUDGET} \
    --no_subprocess

touch ${DONE_TOKEN}
echo "BoltzGen completed at: \$(date)"
SLURM_SCRIPT_EOF

echo "[$(date '+%Y-%m-%d %H:%M:%S')] SLURM script written to $SLURM_SCRIPT"
echo "--- SLURM script contents ---"
cat "$SLURM_SCRIPT"
echo "-----------------------------"

JOB_ID=$(sbatch --parsable "$SLURM_SCRIPT")

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitted SLURM job ${JOB_ID}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Waiting for job to complete..."

# Poll until done or failed
while true; do
    sleep 60
    JOB_STATE=$(squeue -j "$JOB_ID" -h -o "%T" 2>/dev/null || echo "UNKNOWN")

    if [[ -z "$JOB_STATE" || "$JOB_STATE" == "UNKNOWN" ]]; then
        EXIT_CODE=$(sacct -j "$JOB_ID" --format=ExitCode --noheader 2>/dev/null | head -1 | cut -d: -f1 | tr -d ' ')
        if [[ "$EXIT_CODE" == "0" ]] || [[ -f "$DONE_TOKEN" ]]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] BoltzGen job ${JOB_ID} completed successfully"
            break
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] BoltzGen job ${JOB_ID} failed (exit code: ${EXIT_CODE})"
            echo "Check logs: ${BOLTZGEN_DIR}/slurm_${JOB_ID}.out"
            exit 1
        fi
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Job ${JOB_ID} state: ${JOB_STATE}"
done

if [[ ! -f "$DONE_TOKEN" ]]; then
    echo "Error: Done token not found after job completion"
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 1 complete. Output: ${BOLTZGEN_DIR}/output"
