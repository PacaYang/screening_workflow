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
build_binding_yaml() {
    local hotspots="$1"
    python3 <<PYTHON
hotspots = "$hotspots"
from collections import defaultdict
chains = defaultdict(list)
for token in hotspots.split(','):
    token = token.strip()
    if ':' in token:
        chain_id, resnum = token.split(':', 1)
        chains[chain_id.strip()].append(resnum.strip())
    else:
        chains['A'].append(token)

print("      binding_types:")
for chain_id, resnums in sorted(chains.items()):
    print(f"        - chain:")
    print(f"            id: {chain_id}")
    print(f"            binding: {','.join(resnums)}")
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
        build_binding_yaml "$HOTSPOT_RESIDUES"
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
