#!/usr/bin/env bash
# Run binding site analysis for IL6 project proteins via SLURM
# Usage: bash scripts/pipeline/run_binding_site_analysis.sh --task-root ~/Projects/B3/IL6

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

TASK_ROOT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --task-root) TASK_ROOT="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$TASK_ROOT" ]]; then
    echo "Usage: $0 --task-root <path>"
    exit 1
fi

TASK_ROOT="${TASK_ROOT/#\~/$HOME}"

CONDA_BASE="$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")"
CONDA_ENV="boltz_test"
ANALYZE_DIR="${REPO_ROOT}/scripts/scoring"

# Jobs: protein, method, summary_csv_path
declare -A SUMMARY_CSVS=(
    ["IL6RA:AF3"]="${TASK_ROOT}/IL6RA/fine_screening/AF3/summary.csv"
    ["IL6RA:Boltz2"]="${TASK_ROOT}/IL6RA/fine_screening/Boltz2/summary.csv"
    ["IL6RB:AF3"]="${TASK_ROOT}/IL6RB/fine_screening/AF3/summary.csv"
    ["IL6RB:Boltz2"]="${TASK_ROOT}/IL6RB/fine_screening/Boltz2/summary.csv"
    ["IL6RB:RoseTTAFold"]="${TASK_ROOT}/IL6RB/fine_screening/RoseTTAFold/summary.csv"
)

declare -A METHOD_SCRIPTS=(
    ["AF3"]="analyze_af3_binding_sites.py"
    ["Boltz2"]="analyze_boltz2_binding_sites.py"
    ["RoseTTAFold"]="analyze_rosettafold_binding_sites.py"
)

submitted=0
skipped=0

for key in "${!SUMMARY_CSVS[@]}"; do
    PROTEIN="${key%%:*}"
    METHOD="${key##*:}"
    SUMMARY_CSV="${SUMMARY_CSVS[$key]}"
    SCRIPT="${ANALYZE_DIR}/${METHOD_SCRIPTS[$METHOD]}"
    OUTPUT_DIR="${TASK_ROOT}/${PROTEIN}/fine_screening/${METHOD}/binding_site_analysis"
    LOG_DIR="${OUTPUT_DIR}/logs"

    if [[ ! -f "$SUMMARY_CSV" ]]; then
        echo "[SKIP] ${PROTEIN} x ${METHOD}: summary.csv not found at ${SUMMARY_CSV}"
        ((skipped++)) || true
        continue
    fi

    # Check summary.csv has binding_site_center column
    if ! head -1 "$SUMMARY_CSV" | grep -q "binding_site_center"; then
        echo "[SKIP] ${PROTEIN} x ${METHOD}: summary.csv missing binding_site_center column"
        ((skipped++)) || true
        continue
    fi

    mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

    JOB_NAME="bsa_${PROTEIN}_${METHOD}"

    sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --partition=g24
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=${LOG_DIR}/%j.out
#SBATCH --error=${LOG_DIR}/%j.err

set -euo pipefail

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

# Install seaborn if not present
python -c "import seaborn" 2>/dev/null || pip install seaborn -q

python "${SCRIPT}" \
    --summary-csv "${SUMMARY_CSV}" \
    --protein-name "${PROTEIN}" \
    --output-dir "${OUTPUT_DIR}"

echo "Done: ${PROTEIN} x ${METHOD}"
EOF

    echo "[SUBMITTED] ${PROTEIN} x ${METHOD} -> ${OUTPUT_DIR}"
    ((submitted++)) || true
done

echo ""
echo "Submitted: ${submitted} jobs, Skipped: ${skipped} combos"
echo "Monitor with: squeue -u \$USER"
