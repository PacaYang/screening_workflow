#!/bin/bash
# Step 5: Stream AF3 results as jobs complete

set -euo pipefail

# Parse arguments
OUTPUT_DIR=""
PROTEIN_NAME=""
POLL_INTERVAL=60
MAX_ITERATIONS=1000

while [[ $# -gt 0 ]]; do
    case $1 in
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --protein-name) PROTEIN_NAME="$2"; shift 2 ;;
        --poll-interval) POLL_INTERVAL="$2"; shift 2 ;;
        --max-iterations) MAX_ITERATIONS="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$OUTPUT_DIR" || -z "$PROTEIN_NAME" ]]; then
    echo "Error: Missing required arguments"
    exit 1
fi

# Setup paths
AF3_DIR="${OUTPUT_DIR}/04_af3_output"
PRED_DIR="${AF3_DIR}/predictions"
TOKEN_DIR="${PRED_DIR}/token"
STATE_DIR="${AF3_DIR}/.collection_state"
RESULTS_DIR="${OUTPUT_DIR}/05_results"
STREAM_DIR="${RESULTS_DIR}/streaming_updates"

mkdir -p "$STATE_DIR" "$STREAM_DIR"

STATE_FILE="${STATE_DIR}/AF3_peptide.json"

# Initialize state
if [[ ! -f "$STATE_FILE" ]]; then
    echo '{"processed": [], "last_batch": 0}' > "$STATE_FILE"
fi

SCRIPT_DIR="/home/yangl_pacagen_com/screening_workflow/scripts"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting streaming collection"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Polling every ${POLL_INTERVAL}s"

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate base

ITERATION=0
# Total predictions expected = number of AF3 input JSONs. Collection is only
# truly "done" when every prediction has produced a token, regardless of what
# squeue reports (avoids a startup race where jobs haven't registered yet).
EXPECTED_TOTAL=$(find "${OUTPUT_DIR}/03_af3_input" -name "*.json" 2>/dev/null | wc -l)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Expecting ${EXPECTED_TOTAL} predictions"

while [[ $ITERATION -lt $MAX_ITERATIONS ]]; do
    ITERATION=$((ITERATION + 1))

    # Find new completed predictions
    NEW_TOKENS=($(comm -13 \
        <(jq -r '.processed[]' "$STATE_FILE" | sort) \
        <(find "$TOKEN_DIR" -name "*.token" -exec basename {} .token \; | sort)))

    if [[ ${#NEW_TOKENS[@]} -eq 0 ]]; then
        PROCESSED_COUNT=$(jq -r '.processed | length' "$STATE_FILE")

        # Done only when all expected predictions are collected.
        if [[ "$EXPECTED_TOTAL" -gt 0 && "$PROCESSED_COUNT" -ge "$EXPECTED_TOTAL" ]]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] All ${EXPECTED_TOTAL} predictions collected"
            break
        fi

        # Match our AF3 jobs by name prefix (squeue -n needs exact names, so
        # filter the full list instead of relying on a glob).
        RUNNING_JOBS=$(squeue -u "$USER" --noheader -o "%j" 2>/dev/null | grep -c "^${PROTEIN_NAME}_af3_")

        if [[ $RUNNING_JOBS -eq 0 ]]; then
            # No new tokens AND no running jobs. If nothing has been collected
            # yet, the jobs may not have registered — wait a few cycles before
            # giving up; otherwise the run genuinely produced no more results.
            if [[ "$PROCESSED_COUNT" -eq 0 && "$ITERATION" -lt 5 ]]; then
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] No jobs/results yet (startup), waiting..."
                sleep $POLL_INTERVAL
                continue
            fi
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] No running jobs; collected ${PROCESSED_COUNT}/${EXPECTED_TOTAL}"
            break
        fi

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Waiting for results... ($RUNNING_JOBS jobs running, ${PROCESSED_COUNT}/${EXPECTED_TOTAL} done)"
        sleep $POLL_INTERVAL
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found ${#NEW_TOKENS[@]} new results"

    # Process new results
    BATCH_NUM=$(jq -r '.last_batch' "$STATE_FILE")
    BATCH_NUM=$((BATCH_NUM + 1))
    BATCH_FILE="${STREAM_DIR}/batch_${BATCH_NUM}.csv"

    echo "peptide_id,peptide_sequence,backbone_id,pTM,ipTM,pLDDT,num_contacts,ranking_score" > "$BATCH_FILE"

    for token_name in "${NEW_TOKENS[@]}"; do
        # Find the corresponding CIF file. AF3 writes the top-ranked model as
        # "<name>_model.cif" at the top of the prediction dir (per-seed samples
        # live in seed-*/ subdirs as "<name>_seed-*_sample-*_model.cif").
        CIF_FILE="${PRED_DIR}/${token_name}/${token_name}_model.cif"
        if [[ ! -f "$CIF_FILE" ]]; then
            # Fallback: first ranked model anywhere under the prediction dir
            CIF_FILE=$(find "$PRED_DIR/${token_name}" -maxdepth 1 -name "*_model.cif" 2>/dev/null | head -1)
        fi

        if [[ -z "$CIF_FILE" || ! -f "$CIF_FILE" ]]; then
            echo "Warning: CIF not found for $token_name"
            continue
        fi

        # Parse peptide info from name (format: PROTEIN_backboneID_seqXXX)
        BACKBONE_ID=$(echo "$token_name" | sed -E "s/${PROTEIN_NAME}_(.*)_seq[0-9]+/\1/")

        # Extract scores using Python
        python3 <<PYTHON_SCRIPT
import json
import sys
from pathlib import Path

try:
    from Bio.PDB import MMCIFParser
    from Bio.Data.PDBData import protein_letters_3to1_extended as THREE_TO_ONE
    import numpy as np

    cif_file = '${CIF_FILE}'
    pred_dir = Path('${PRED_DIR}/${token_name}')

    # pTM / ipTM come from AF3's summary_confidences.json, not the CIF.
    ptm, iptm = 0.0, 0.0
    summary = pred_dir / '${token_name}_summary_confidences.json'
    if summary.is_file():
        sc = json.load(open(summary))
        ptm = sc.get('ptm') or 0.0
        iptm = sc.get('iptm') or 0.0

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('peptide', cif_file)
    model = structure[0]

    # The peptide is the SHORTEST chain (target protomers are longer).
    chain_lengths = {ch.id: sum(1 for r in ch if r.id[0] == ' ') for ch in model}
    pep_chain_id = min(chain_lengths, key=chain_lengths.get)
    target_chain_ids = [c for c in chain_lengths if c != pep_chain_id]

    # pLDDT from B-factors (AF3 stores per-atom pLDDT in the B-factor column).
    plddt_values = [atom.bfactor for ch in model for res in ch for atom in res]
    plddt = float(np.mean(plddt_values)) if plddt_values else 0.0

    # Count peptide residues contacting ANY target chain (min heavy-atom dist <= 5A)
    pep_chain = model[pep_chain_id]
    target_atoms = [atom.coord for tcid in target_chain_ids for res in model[tcid]
                    if res.id[0] == ' ' for atom in res]
    target_atoms = np.array(target_atoms) if target_atoms else np.empty((0, 3))

    num_contacts = 0
    pep_seq = ''
    for res in pep_chain:
        if res.id[0] != ' ':
            continue
        try:
            pep_seq += THREE_TO_ONE.get(res.get_resname(), 'X')
        except Exception:
            pep_seq += 'X'
        if len(target_atoms):
            for atom in res:
                d = np.linalg.norm(target_atoms - atom.coord, axis=1)
                if d.min() <= 5.0:
                    num_contacts += 1
                    break

    ranking_score = 0.4 * iptm + 0.3 * ptm + 0.3 * (plddt / 100.0)

    print(f"${token_name},{pep_seq},${BACKBONE_ID},{ptm:.4f},{iptm:.4f},{plddt:.2f},{num_contacts},{ranking_score:.4f}")

except Exception as e:
    print(f"Error processing ${token_name}: {e}", file=sys.stderr)
PYTHON_SCRIPT

    done >> "$BATCH_FILE" 2>/dev/null

    # Update state
    jq --arg tokens "$(printf '%s\n' "${NEW_TOKENS[@]}" | jq -R . | jq -s .)" \
       --arg batch "$BATCH_NUM" \
       '.processed += ($tokens | fromjson) | .last_batch = ($batch | tonumber)' \
       "$STATE_FILE" > "${STATE_FILE}.tmp" && mv "${STATE_FILE}.tmp" "$STATE_FILE"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Batch ${BATCH_NUM} saved to ${BATCH_FILE}"

    sleep $POLL_INTERVAL
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Streaming collection complete"
N_PROCESSED=$(jq -r '.processed | length' "$STATE_FILE")
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Total processed: ${N_PROCESSED}"
