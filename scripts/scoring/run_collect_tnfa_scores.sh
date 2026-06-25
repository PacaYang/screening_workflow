#!/bin/bash
#
# Collect fine-screening scores for TNFa (trimer) across all 4 methods.
# Extracts each method's tarballs, runs the per-method scorers, and merges
# into combined_summary.csv. Idempotent: skips a method if its summary exists.
#
set -e

source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

TASK_ROOT="/shared/B3/TNFa_Jun10"
PROTEINS="TNFa"
SCRIPTS_ROOT="/home/ubuntu/screening_workflow/scripts"

GENERAL_PYTHON="/home/ubuntu/miniconda3/envs/general/bin/python"
RFAA_PYTHON="/home/ubuntu/miniconda3/envs/RFAA/bin/python"

log_info() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"; }

for protein in $PROTEINS; do
    base="${TASK_ROOT}/${protein}"
    selected="${base}/initial_screening/selected.csv"
    log_info "=== Processing ${protein} ==="

    # ---- AF3 ----
    af3_dir="${base}/fine_screening/AF3"
    af3_out="${af3_dir}/output"
    if [ ! -f "${af3_dir}/summary.csv" ] && [ -d "${af3_out}" ]; then
        log_info "AF3: extracting tarballs..."
        for f in "${af3_out}"/batch_*.tar.gz; do
            [ -f "$f" ] && tar -xzf "$f" -C "${af3_out}/" 2>/dev/null || true
        done
        log_info "AF3: running scorer..."
        $GENERAL_PYTHON "${SCRIPTS_ROOT}/scoring/af3_scores.py" \
            --af3-results-folder "${af3_out}" \
            --output-dir "${af3_dir}" || log_info "AF3 scorer failed for ${protein}"
    else
        log_info "AF3: summary exists or no output, skipping."
    fi

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
            --smiles "${selected}" \
            --smiles-col "SMILES" || log_info "Boltz2 scorer failed for ${protein}"
    else
        log_info "Boltz2: summary exists or no output, skipping."
    fi

    # ---- RoseTTAFold ----
    # Guard: don't score a partial RFAA run. Skip if RFAA jobs are still queued,
    # since the per-method summary.csv is written once and then treated as final.
    rfaa_dir="${base}/fine_screening/RoseTTAFold"
    rfaa_out="${rfaa_dir}/protein_ligand/output"
    RFAA_QUEUED=$(squeue -u "$(whoami)" -h -o "%j" 2>/dev/null | grep -ic -E "rfaa|rosetta" || true)
    if [ "${RFAA_QUEUED:-0}" -gt 0 ]; then
        log_info "RFAA: ${RFAA_QUEUED} jobs still in queue, deferring RFAA scoring."
    elif [ ! -f "${rfaa_dir}/summary.csv" ] && [ -d "${rfaa_out}" ]; then
        log_info "RFAA: extracting tarballs..."
        for f in "${rfaa_out}"/batch_*.tar.gz; do
            [ -f "$f" ] && tar -xzf "$f" -C "${rfaa_out}/" 2>/dev/null || true
        done
        log_info "RFAA: running scorer..."
        $RFAA_PYTHON "${SCRIPTS_ROOT}/scoring/rosettafold_scores.py" \
            --rfaa-results-folder "${rfaa_out}" \
            --protein-name "${protein}" \
            --output-dir "${rfaa_dir}" || log_info "RFAA scorer failed for ${protein}"
    else
        log_info "RFAA: summary exists or no output, skipping."
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
    else
        log_info "Vina: results exist or no output, skipping."
    fi

    # ---- Merge ----
    log_info "Merging ${protein}..."
    $GENERAL_PYTHON - "$base" <<'PYEOF'
import pandas as pd, os, sys
base = sys.argv[1]
out_path = os.path.join(base, "fine_screening", "combined_summary.csv")
def load(path, prefix):
    if not os.path.exists(path):
        print(f"  MISSING: {path}"); return None
    try: df = pd.read_csv(path)
    except Exception as e:
        print(f"  SKIP (unreadable): {path} — {e}"); return None
    if df.empty or len(df.columns) == 0:
        print(f"  SKIP (empty): {path}"); return None
    df = df.loc[:, ~df.columns.str.match(r'^Unnamed')]
    smi = next((c for c in df.columns if c.upper() == 'SMILES'), None)
    if smi is None:
        print(f"  No SMILES col in {path}: {list(df.columns)}"); return None
    df = df.rename(columns={smi: 'SMILES'})
    df = df.rename(columns={c: f"{prefix}_{c}" for c in df.columns if c != 'SMILES'})
    print(f"  Loaded {len(df)} rows from {os.path.basename(os.path.dirname(path))}")
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
    print(f"  Combined summary ({len(combined)} rows) -> {out_path}")
PYEOF
    log_info "Done with ${protein}"
done
log_info "All complete."
