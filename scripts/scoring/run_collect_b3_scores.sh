#!/bin/bash
#SBATCH --job-name=collect_b3_scores
#SBATCH --output=/shared/B3/IL6_Apr13/collect_scores_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

set -e

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="/shared/B3/IL6_Apr13"
PROTEINS="IL6RA IL6RB"
SCRIPTS_ROOT="/home/ubuntu/screening_workflow/scripts"

GENERAL_PYTHON="/home/ubuntu/miniconda3/envs/general/bin/python"
RFAA_PYTHON="/home/ubuntu/miniconda3/envs/RFAA/bin/python"

log_info() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; }

for protein in $PROTEINS; do
    base="${TASK_ROOT}/${protein}"
    log_info "=== Processing ${protein} ==="

    # ---- Boltz2 ----
    boltz_out="${base}/fine_screening/Boltz2/output"
    boltz_dir="${base}/fine_screening/Boltz2"
    selected_af3="${base}/fine_screening/AF3_first/selected_af3.csv"

    if [ ! -f "${boltz_dir}/summary.csv" ]; then
        log_info "Extracting Boltz2 tar.gz archives for ${protein}..."
        for f in "${boltz_out}"/batch_*.tar.gz; do
            [ -f "$f" ] && tar -xzf "$f" -C "${boltz_out}/"
        done

        log_info "Running boltz2_scores.py for ${protein}..."
        $GENERAL_PYTHON "${SCRIPTS_ROOT}/scoring/boltz2_scores.py" \
            --boltz-results-folder "${boltz_out}" \
            --output-dir "${boltz_dir}" \
            --smiles "${selected_af3}" \
            --smiles-col "SMILES"
        log_info "Boltz2 summary written to ${boltz_dir}/summary.csv"
    else
        log_info "Boltz2 summary already exists for ${protein}, skipping."
    fi

    # ---- RoseTTAFold (needs torch, use RFAA env) ----
    rfaa_dir="${base}/fine_screening/RoseTTAFold"
    rfaa_out="${rfaa_dir}/protein_ligand/output"

    if [ ! -f "${rfaa_dir}/summary.csv" ]; then
        log_info "Pre-extracting RoseTTAFold archives for ${protein} (faster than Python tarfile)..."
        for f in "${rfaa_out}"/batch_*.tar.gz; do
            [ -f "$f" ] && tar -xzf "$f" -C "${rfaa_out}/"
        done
        log_info "Running rosettafold_scores.py for ${protein}..."
        $RFAA_PYTHON "${SCRIPTS_ROOT}/scoring/rosettafold_scores.py" \
            --rfaa-results-folder "${rfaa_out}" \
            --protein-name "${protein}" \
            --output-dir "${rfaa_dir}"
        log_info "RoseTTAFold summary written to ${rfaa_dir}/summary.csv"
    else
        log_info "RoseTTAFold summary already exists for ${protein}, skipping."
    fi

    # ---- Vina ----
    vina_dir="${base}/fine_screening/Vina"
    vina_out="${vina_dir}/output"
    vina_input="${vina_dir}/input"

    if [ ! -f "${vina_dir}/results.csv" ]; then
        log_info "Running vina_scores.py for ${protein}..."
        $GENERAL_PYTHON "${SCRIPTS_ROOT}/scoring/vina_scores.py" \
            --vina-results-folder "${vina_out}" \
            --output-dir "${vina_dir}" \
            --input-dir "${vina_input}"
        log_info "Vina results written to ${vina_dir}/results.csv"
    else
        log_info "Vina results already exist for ${protein}, skipping."
    fi

    # ---- Merge all methods into combined_summary.csv ----
    log_info "Merging all method scores for ${protein}..."
    $GENERAL_PYTHON - <<PYEOF
import pandas as pd
import os

base = "${base}"
out_path = os.path.join(base, "fine_screening", "combined_summary.csv")

def load(path, prefix):
    if not os.path.exists(path):
        print(f"  MISSING: {path}")
        return None
    try:
        df = pd.read_csv(path)
    except Exception as e:
        print(f"  SKIP (unreadable): {path} — {e}")
        return None
    if df.empty or len(df.columns) == 0:
        print(f"  SKIP (empty): {path}")
        return None
    # Drop unnamed index columns
    df = df.loc[:, ~df.columns.str.match(r'^Unnamed')]
    # Identify SMILES column (normalize name)
    smiles_col = next((c for c in df.columns if c.upper() == 'SMILES'), None)
    if smiles_col is None:
        print(f"  No SMILES column in {path}, columns: {list(df.columns)}")
        return None
    df = df.rename(columns={smiles_col: 'SMILES'})
    # Prefix non-SMILES columns
    rename = {c: f"{prefix}_{c}" for c in df.columns if c != 'SMILES'}
    df = df.rename(columns=rename)
    print(f"  Loaded {len(df)} rows from {path}")
    return df

af3   = load(os.path.join(base, "fine_screening", "AF3", "summary.csv"), "af3")
boltz = load(os.path.join(base, "fine_screening", "Boltz2", "summary.csv"), "boltz2")
rfaa  = load(os.path.join(base, "fine_screening", "RoseTTAFold", "summary.csv"), "rfaa")
vina  = load(os.path.join(base, "fine_screening", "Vina", "results.csv"), "vina")

dfs = [d for d in [af3, boltz, rfaa, vina] if d is not None]
if not dfs:
    print("No data found, skipping merge.")
else:
    combined = dfs[0]
    for d in dfs[1:]:
        combined = combined.merge(d, on='SMILES', how='outer')
    combined.to_csv(out_path, index=False)
    print(f"  Combined summary ({len(combined)} rows) written to {out_path}")
PYEOF

    log_info "Done with ${protein}"
done

log_info "All proteins complete."
