#!/bin/bash
#SBATCH --job-name=collect_woundcare_scores
#SBATCH --output=/shared/cuteness_woundcare/logs/collect_scores_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

set -e

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="/shared/cuteness_woundcare"
SCRIPTS_ROOT="/home/ubuntu/screening_workflow/scripts"

GENERAL_PYTHON="/home/ubuntu/miniconda3/envs/general/bin/python"
RFAA_PYTHON="/home/ubuntu/miniconda3/envs/RFAA/bin/python"

EXCLUDE_REGEX='^(JAK1|JAK2|JAK3|TYK2|JAK_input|Input|logs)$'
RFAA_SKIP=" RPTOR PIK3CB MLXIPL FASN ACACA LATS1 PDE3B PLOD2 PLOD1 PDE4D MTOR RICTOR ADCY1 ADCY3 ADCY5 ADCY6 ADCY9 ITGA1 ITGA2 SREBF1 "

mkdir -p "${TASK_ROOT}/logs"

log_info() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"; }

PROTEINS=()
for entry in "${TASK_ROOT}"/*; do
    [ -d "$entry" ] || continue
    name=$(basename "$entry")
    if [[ "$name" =~ $EXCLUDE_REGEX ]]; then continue; fi
    if [ -d "${entry}/fine_screening" ]; then PROTEINS+=("$name"); fi
done

log_info "Found ${#PROTEINS[@]} proteins"

for protein in "${PROTEINS[@]}"; do
    base="${TASK_ROOT}/${protein}"
    log_info "=== ${protein} ==="

    selected_smiles="${base}/initial_screening/selected_all21.csv"

    # ---- Boltz2 ----
    boltz_dir="${base}/fine_screening/Boltz2"
    boltz_out="${boltz_dir}/output"
    if [ ! -f "${boltz_dir}/summary.csv" ] && [ -d "${boltz_out}" ]; then
        log_info "Boltz2: extracting tarballs..."
        for f in "${boltz_out}"/batch_*.tar.gz; do
            [ -f "$f" ] && tar -xzf "$f" -C "${boltz_out}/" 2>/dev/null || true
        done
        log_info "Boltz2: running scorer..."
        $GENERAL_PYTHON "${SCRIPTS_ROOT}/scoring/boltz2_scores.py" \
            --boltz-results-folder "${boltz_out}" \
            --output-dir "${boltz_dir}" \
            --smiles "${selected_smiles}" \
            --smiles-col "SMILES" || log_info "Boltz2 scorer failed for ${protein}"
    fi

    # ---- RoseTTAFold (skip prefold-failed proteins) ----
    if [[ ! "$RFAA_SKIP" =~ " ${protein} " ]]; then
        rfaa_dir="${base}/fine_screening/RoseTTAFold"
        rfaa_out="${rfaa_dir}/protein_ligand/output"
        if [ ! -f "${rfaa_dir}/summary.csv" ] && [ -d "${rfaa_out}" ]; then
            log_info "RFAA: extracting tarballs..."
            for f in "${rfaa_out}"/batch_*.tar.gz; do
                [ -f "$f" ] && tar -xzf "$f" -C "${rfaa_out}/" 2>/dev/null || true
            done
            log_info "RFAA: running scorer..."
            $RFAA_PYTHON "${SCRIPTS_ROOT}/scoring/rosettafold_scores.py" \
                --rfaa-results-folder "${rfaa_out}" \
                --protein-name "${protein}" \
                --output-dir "${rfaa_dir}" || log_info "RFAA scorer failed for ${protein}"
        fi
    fi

    # ---- Vina ----
    vina_dir="${base}/fine_screening/Vina"
    vina_out="${vina_dir}/output"
    vina_input="${vina_dir}/input"
    if [ ! -f "${vina_dir}/results.csv" ] && [ -d "${vina_out}" ]; then
        log_info "Vina: running scorer..."
        $GENERAL_PYTHON "${SCRIPTS_ROOT}/scoring/vina_scores.py" \
            --vina-results-folder "${vina_out}" \
            --output-dir "${vina_dir}" \
            --input-dir "${vina_input}" || log_info "Vina scorer failed for ${protein}"
    fi

    # ---- Merge ----
    log_info "Merging ${protein}..."
    $GENERAL_PYTHON - <<PYEOF
import pandas as pd, os
base = "${base}"
out_path = os.path.join(base, "fine_screening", "combined_summary.csv")
def load(path, prefix):
    if not os.path.exists(path): return None
    try: df = pd.read_csv(path)
    except: return None
    if df.empty or len(df.columns) == 0: return None
    df = df.loc[:, ~df.columns.str.match(r'^Unnamed')]
    smi = next((c for c in df.columns if c.upper() == 'SMILES'), None)
    if smi is None: return None
    df = df.rename(columns={smi: 'SMILES'})
    df = df.rename(columns={c: f"{prefix}_{c}" for c in df.columns if c != 'SMILES'})
    return df
af3   = load(os.path.join(base, "fine_screening", "AF3", "summary.csv"), "af3")
boltz = load(os.path.join(base, "fine_screening", "Boltz2", "summary.csv"), "boltz2")
rfaa  = load(os.path.join(base, "fine_screening", "RoseTTAFold", "summary.csv"), "rfaa")
vina  = load(os.path.join(base, "fine_screening", "Vina", "results.csv"), "vina")
dfs = [d for d in [af3, boltz, rfaa, vina] if d is not None]
if dfs:
    combined = dfs[0]
    for d in dfs[1:]:
        combined = combined.merge(d, on='SMILES', how='outer')
    combined.to_csv(out_path, index=False)
    print(f"  combined ({len(combined)} rows): {out_path}")
PYEOF
done

log_info "All proteins complete."
