#!/bin/bash
# Design 30-AA peptide binders against IL13RA1 (PDB 3BPN, chain C).
#
# Target epitope (chain C of 3BPN, IL-13 receptor alpha-1 extracellular domain):
#   R256, C257, K318, L319, C320, Y321, E322
#
# Two design engines, each producing 900 peptides (1800 total), parallelised
# across the cluster as SLURM job arrays -- "one design unit per node":
#
#   RFDiffusion : array of 30 GPU tasks, each = 1 backbone x 30 ProteinMPNN seqs.
#                 Each task runs the FULL run_peptide_pipeline.sh end-to-end
#                 (RFDiffusion -> ProteinMPNN -> AF3 -> rank) in its own dir.
#                 30 tasks x 30 = 900 peptides.
#   BoltzGen    : array of 30 CPU tasks, each designs 30 peptides. BoltzGen's
#                 own GPU work (design + AF3) self-submits as sub-jobs, so these
#                 orchestrator tasks need no GPU. 30 tasks x 30 = 900 peptides.
#
# Every task is fully independent: its own output subdir, its own .pipeline_state
# and its own ranked_peptides.csv. Merge afterwards with the printed find/sort.
#
# Both cysteines in the epitope (C257, C320) are disulfide-bonded in the native
# receptor fold; they are used only as CA hotspots here, not covalent anchors.
#
# NOTE: RFDiffusion array tasks hold their GPU while their AF3 sub-jobs run and
# are polled (steps 4-5). That is the cost of "fully independent per backbone".
#
# Usage:
#   bash design_il13ra1_3bpn.sh                       # submit both engine arrays
#   bash design_il13ra1_3bpn.sh --engine rfdiffusion  # only RFDiffusion -> 900
#   bash design_il13ra1_3bpn.sh --engine boltzgen     # only BoltzGen    -> 900
#   bash design_il13ra1_3bpn.sh --partition g232      # GPU partition for RFD array
#   bash design_il13ra1_3bpn.sh --max-parallel 15     # throttle concurrent tasks
#   bash design_il13ra1_3bpn.sh --dry-run             # write scripts, don't submit

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Target definition (verified against 3BPN.pdb: chain C = IL13RA1) ---------
PROTEIN="IL13RA1"
PDB_ID="3BPN"
PEPTIDE_LENGTH=30
# Chains in 3BPN: A = IL-4, B = IL-4Ra, C = IL13RA1 (our target).
HOTSPOT_RESIDUES="C:256,C:257,C:318,C:319,C:320,C:321,C:322"

# --- Parallel layout: 30 tasks per engine, 30 peptides each -> 900 each --------
RFD_N_TASKS=30        # RFDiffusion array size (one backbone per task)
RFD_SEQS_PER_BB=30    # ProteinMPNN sequences per backbone -> 1*30 = 30 / task
BOLTZ_N_TASKS=30      # BoltzGen array size
BOLTZ_PER_TASK=30     # BoltzGen designs per task          ->         30 / task
AF3_N_BATCHES=5       # AF3 SLURM batches per task (only 30 peptides per task)

# --- Paths --------------------------------------------------------------------
PROTEIN_PDB="${IL13RA1_PDB:-/home/yangl_pacagen_com/Projects/B3/IL4_peptide/3BPN.pdb}"
BASE_OUTDIR="${IL13RA1_OUTDIR:-/home/yangl_pacagen_com/Projects/B3/IL4_peptide/il13ra1_peptides}"

# --- SLURM defaults -----------------------------------------------------------
# Per-node limits: g24 = 2 CPU/15G, g212 = 6 CPU/47G, g232 = 16 CPU/124G (all gpu:1).
PARTITION="g232"      # GPU partition for RFDiffusion array (needs the big nodes)
CPU_PARTITION="g24"   # BoltzGen CPU orchestrator array (tiny: 2 CPU/8G fits g24)
GPU_TIME="48:00:00"
GPU_MEM="64G"
GPU_CPUS=8
CPU_TIME="48:00:00"
CPU_MEM="8G"
MAX_PARALLEL=0        # 0 = no %N throttle on the arrays
CONDA_SH="/home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh"

# --- Options ------------------------------------------------------------------
ENGINE="both"
DRY_RUN=0
PASSTHRU=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --engine)        ENGINE="$2"; shift 2 ;;
        --partition)     PARTITION="$2"; shift 2 ;;
        --cpu-partition) CPU_PARTITION="$2"; shift 2 ;;
        --max-parallel)  MAX_PARALLEL="$2"; shift 2 ;;
        --dry-run)       DRY_RUN=1; shift ;;
        *) PASSTHRU+=("$1"); shift ;;   # forwarded to run_*_pipeline.sh
    esac
done

case "$ENGINE" in
    both|rfdiffusion|boltzgen) ;;
    *) echo "Error: --engine must be one of: both, rfdiffusion, boltzgen"; exit 1 ;;
esac

# Array %N throttle suffix, if requested.
ARRAY_THROTTLE=""
[[ "$MAX_PARALLEL" -gt 0 ]] && ARRAY_THROTTLE="%${MAX_PARALLEL}"

# Fetch the structure if it is not already present.
if [[ ! -f "$PROTEIN_PDB" ]]; then
    echo "[design] ${PDB_ID}.pdb not found at ${PROTEIN_PDB}; downloading from RCSB..."
    mkdir -p "$(dirname "$PROTEIN_PDB")"
    curl -sfL "https://files.rcsb.org/download/${PDB_ID}.pdb" -o "$PROTEIN_PDB"
fi

echo "[design] Target:          ${PROTEIN} (${PDB_ID} chain C)"
echo "[design] Peptide length:  ${PEPTIDE_LENGTH} AA"
echo "[design] Hotspot epitope: ${HOTSPOT_RESIDUES}"
echo "[design] PDB:             ${PROTEIN_PDB}"
echo "[design] Engine(s):       ${ENGINE}"
echo "[design] Base output:     ${BASE_OUTDIR}"
[[ $DRY_RUN -eq 1 ]] && echo "[design] DRY RUN: scripts written but not submitted"

PASSTHRU_STR=""
[[ ${#PASSTHRU[@]} -gt 0 ]] && PASSTHRU_STR="$(printf '%q ' "${PASSTHRU[@]}")"

submit_array() {   # $1 = path to the array script (already written)
    local script=$1
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "[design] (dry-run) would submit: sbatch ${script}"
        return
    fi
    local jid; jid=$(sbatch --parsable "$script")
    echo "[design] Submitted array job ${jid}: ${script}"
}

# ---------------------------------------------------------------------------
# RFDiffusion: 30 GPU tasks, each a full pipeline for one backbone.
# ---------------------------------------------------------------------------
build_rfdiffusion_array() {
    local root="${BASE_OUTDIR}/rfdiffusion"
    mkdir -p "${root}/logs"
    local script="${root}/rfd_array.slurm"
    cat > "$script" <<SLURM
#!/bin/bash
#SBATCH --job-name=${PROTEIN}_rfd
#SBATCH --output=${root}/logs/task_%a_%A.out
#SBATCH --error=${root}/logs/task_%a_%A.err
#SBATCH --array=0-$((RFD_N_TASKS-1))${ARRAY_THROTTLE}
#SBATCH --time=${GPU_TIME}
#SBATCH --mem=${GPU_MEM}
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=${GPU_CPUS}
#SBATCH --partition=${PARTITION}

set -euo pipefail
TID=\$(printf '%02d' "\$SLURM_ARRAY_TASK_ID")
TASK_OUT="${root}/bb_\${TID}"
echo "RFD task \$SLURM_ARRAY_TASK_ID on \$(hostname) -> \${TASK_OUT} at \$(date)"

source ${CONDA_SH}
bash "${SCRIPT_DIR}/run_peptide_pipeline.sh" \\
    --protein "${PROTEIN}" \\
    --protein-pdb "${PROTEIN_PDB}" \\
    --peptide-length ${PEPTIDE_LENGTH} \\
    --output-dir "\${TASK_OUT}" \\
    --hotspot-residues "${HOTSPOT_RESIDUES}" \\
    --n-designs 1 \\
    --n-seqs ${RFD_SEQS_PER_BB} \\
    --n-batches ${AF3_N_BATCHES} \\
    ${PASSTHRU_STR}
echo "RFD task \$SLURM_ARRAY_TASK_ID done at \$(date)"
SLURM
    echo "$script"
}

# ---------------------------------------------------------------------------
# BoltzGen: 30 CPU orchestrator tasks (GPU work self-submits as sub-jobs).
# ---------------------------------------------------------------------------
build_boltzgen_array() {
    local root="${BASE_OUTDIR}/boltzgen"
    mkdir -p "${root}/logs"
    local script="${root}/boltz_array.slurm"
    cat > "$script" <<SLURM
#!/bin/bash
#SBATCH --job-name=${PROTEIN}_boltz
#SBATCH --output=${root}/logs/task_%a_%A.out
#SBATCH --error=${root}/logs/task_%a_%A.err
#SBATCH --array=0-$((BOLTZ_N_TASKS-1))${ARRAY_THROTTLE}
#SBATCH --time=${CPU_TIME}
#SBATCH --mem=${CPU_MEM}
#SBATCH --cpus-per-task=2
#SBATCH --partition=${CPU_PARTITION}

set -euo pipefail
TID=\$(printf '%02d' "\$SLURM_ARRAY_TASK_ID")
TASK_OUT="${root}/chunk_\${TID}"
echo "BoltzGen task \$SLURM_ARRAY_TASK_ID on \$(hostname) -> \${TASK_OUT} at \$(date)"

source ${CONDA_SH}
bash "${SCRIPT_DIR}/boltzgen/run_boltzgen_af3_pipeline.sh" \\
    --protein "${PROTEIN}" \\
    --protein-pdb "${PROTEIN_PDB}" \\
    --peptide-length ${PEPTIDE_LENGTH} \\
    --output-dir "\${TASK_OUT}" \\
    --hotspot-residues "${HOTSPOT_RESIDUES}" \\
    --n-designs ${BOLTZ_PER_TASK} \\
    --n-batches ${AF3_N_BATCHES} \\
    ${PASSTHRU_STR}
echo "BoltzGen task \$SLURM_ARRAY_TASK_ID done at \$(date)"
SLURM
    echo "$script"
}

if [[ "$ENGINE" == "rfdiffusion" || "$ENGINE" == "both" ]]; then
    echo ""
    echo "[design] RFDiffusion: ${RFD_N_TASKS} GPU tasks x ${RFD_SEQS_PER_BB} seqs = $((RFD_N_TASKS*RFD_SEQS_PER_BB)) peptides -> ${PARTITION}"
    submit_array "$(build_rfdiffusion_array)"
fi

if [[ "$ENGINE" == "boltzgen" || "$ENGINE" == "both" ]]; then
    echo ""
    echo "[design] BoltzGen: ${BOLTZ_N_TASKS} CPU tasks x ${BOLTZ_PER_TASK} designs = $((BOLTZ_N_TASKS*BOLTZ_PER_TASK)) peptides (GPU self-submitted)"
    submit_array "$(build_boltzgen_array)"
fi

cat <<INFO

[design] Monitor:  squeue -u \$USER
[design] Per-task results (one CSV each):
           ${BASE_OUTDIR}/rfdiffusion/bb_*/05_results/ranked_peptides.csv
           ${BASE_OUTDIR}/boltzgen/chunk_*/05_results/ranked_peptides.csv
[design] Merge all 1800 into a single ranking (header from the first file):
           find ${BASE_OUTDIR} -path '*/05_results/ranked_peptides.csv' \\
             | xargs awk 'FNR==1 && NR!=1{next} {print}' > ${BASE_OUTDIR}/all_ranked_peptides.csv
INFO
