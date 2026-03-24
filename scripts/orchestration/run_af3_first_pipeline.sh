#!/bin/bash
#
# AF3-first fine-screening pipeline
# 1) Write prefold inputs (AF3/RoseTTAFold/Boltz2) + prepare all-compound selected lists
# 2) Run prefold jobs (AF3/Boltz2/RoseTTAFold) + Vina receptor preparation
# 3) Run AF3 on all compounds
# 4) Collect AF3 scores and select top AF3 compounds
# 5) Run Vina/Boltz2/RoseTTAFold on AF3-selected compounds
# 6) Collect results
#

set -e

TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
PROTEINS="${MASTER_PROTEINS:-JAK1JH1}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SELF_SCRIPT="${SCRIPTS_ROOT}/orchestration/run_af3_first_pipeline.sh"

START_FROM=1
STOP_AFTER=6
POLL_INTERVAL=300
DRY_RUN=0
CURRENT_STAGE=0

COMPOUNDS_CSV="${TASK_ROOT}/Input/compounds_smiles.csv"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"
SMILES_COL="SMILES"
VINA_CONDA_ENV="${VINA_CONDA_ENV:-vina_new}"

SELECTED_REL_ALL="fine_screening/AF3_first/selected_all.csv"
SELECTED_REL_AF3="fine_screening/AF3_first/selected_af3.csv"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

PIPELINE_LOG=""

log_info() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
    echo "=========================================="
    echo "$msg"
    echo "=========================================="
    echo ""
    [ -n "$PIPELINE_LOG" ] && echo "$msg" >> "$PIPELINE_LOG"
}

log_error() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*"
    echo "=========================================="
    echo -e "${RED}${msg}${NC}"
    echo "=========================================="
    echo ""
    [ -n "$PIPELINE_LOG" ] && echo "$msg" >> "$PIPELINE_LOG"
}

check_script_exists() {
    local script=$1
    if [ ! -f "$script" ]; then
        log_error "Script not found: $script"
        return 1
    fi
    return 0
}

cleanup() {
    echo ""
    log_error "Pipeline interrupted at stage ${CURRENT_STAGE}"
    echo "To resume, run:"
    echo "  $SELF_SCRIPT --start-from ${CURRENT_STAGE} --task-root \"$TASK_ROOT\" --proteins \"$PROTEINS\""
    echo ""
    exit 130
}
trap cleanup SIGINT SIGTERM

wait_for_slurm_jobs() {
    local job_file=$1
    local label=$2

    if [ ! -f "$job_file" ] || [ ! -s "$job_file" ]; then
        log_info "No jobs to wait for ($label)"
        return 0
    fi

    local total
    total=$(wc -l < "$job_file")
    log_info "Waiting for ${total} ${label} SLURM jobs..."

    while true; do
        local still_running=0
        while IFS= read -r job_id; do
            job_id=$(echo "$job_id" | tr -d '[:space:]')
            [ -z "$job_id" ] && continue
            if squeue -j "$job_id" &>/dev/null && squeue -j "$job_id" 2>/dev/null | grep -q "$job_id"; then
                still_running=$((still_running + 1))
            fi
        done < "$job_file"

        if [ "$still_running" -eq 0 ]; then
            log_info "All ${label} jobs completed"
            return 0
        fi

        log_info "${label}: ${still_running}/${total} jobs still running/pending. Next check in ${POLL_INTERVAL}s"
        sleep "$POLL_INTERVAL"
    done
}

collect_all_job_files() {
    local dir_pattern=$1
    local tmp_file
    tmp_file=$(mktemp)

    for protein in $PROTEINS; do
        local jf="${TASK_ROOT}/${protein}/${dir_pattern}/submitted_jobs.txt"
        if [ -f "$jf" ]; then
            cat "$jf" >> "$tmp_file"
        fi
    done
    echo "$tmp_file"
}

run_or_dry() {
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would run: $*"
    else
        "$@"
    fi
}

run_stage() {
    local num=$1
    local name=$2
    local func=$3

    if [ "$num" -lt "$START_FROM" ]; then
        echo "Skipping stage ${num} (${name}) — before --start-from ${START_FROM}"
        echo ""
        return 0
    fi
    if [ "$num" -gt "$STOP_AFTER" ]; then
        echo "Skipping stage ${num} (${name}) — after --stop-after ${STOP_AFTER}"
        echo ""
        return 0
    fi

    CURRENT_STAGE=$num
    log_info "=== Stage ${num}: ${name} ==="

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "[DRY-RUN] Would execute stage ${num}: ${name}"
        $func
        echo ""
        return 0
    fi

    $func

    log_info "=== Stage ${num}: ${name} — complete ==="
}

preflight_check_scripts() {
    local required=(
        "$SCRIPTS_ROOT/prepare_input/run_write_prefold_af3.sh"
        "$SCRIPTS_ROOT/prepare_input/run_write_prefold_boltz2.sh"
        "$SCRIPTS_ROOT/structure_prediction/run_prefold_af3.sh"
        "$SCRIPTS_ROOT/structure_prediction/run_prefold_boltz2.sh"
        "$SCRIPTS_ROOT/structure_prediction/run_rosettafold_prefold.sh"
        "$SCRIPTS_ROOT/prepare_input/prepare_vina_receptor.sh"
        "$SCRIPTS_ROOT/prepare_input/run_write_af3_input.sh"
        "$SCRIPTS_ROOT/structure_prediction/run_af3_batch.sh"
        "$SCRIPTS_ROOT/scoring/af3_scores.py"
        "$SCRIPTS_ROOT/scoring/select_af3_top.py"
        "$SCRIPTS_ROOT/prepare_input/run_write_boltz2_input.sh"
        "$SCRIPTS_ROOT/docking/run_vina_batch.sh"
        "$SCRIPTS_ROOT/structure_prediction/run_boltz2_batch.sh"
        "$SCRIPTS_ROOT/structure_prediction/run_rosettafold_batch.sh"
        "$SCRIPTS_ROOT/scoring/vina_scores.py"
        "$SCRIPTS_ROOT/scoring/boltz2_scores.py"
        "$SCRIPTS_ROOT/scoring/rosettafold_scores.py"
    )

    local missing=0
    local script
    for script in "${required[@]}"; do
        if ! check_script_exists "$script"; then
            missing=1
        fi
    done

    if [ ! -f "$COMPOUNDS_CSV" ]; then
        log_error "Compounds CSV not found: $COMPOUNDS_CSV"
        missing=1
    fi
    if [ ! -f "$SEQS_CSV" ]; then
        log_error "Sequences CSV not found: $SEQS_CSV"
        missing=1
    fi

    [ "$missing" -eq 0 ]
}

prepare_all_compounds_for_protein() {
    local protein=$1
    local outdir="${TASK_ROOT}/${protein}/fine_screening/AF3_first"
    local out_csv="${outdir}/selected_all.csv"

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would write AF3-all compounds: ${out_csv}"
        return 0
    fi

    mkdir -p "$outdir"
    python3 - "$COMPOUNDS_CSV" "$SMILES_COL" "$out_csv" <<'PY'
import sys
import pandas as pd

input_csv, smiles_col, output_csv = sys.argv[1:4]
df = pd.read_csv(input_csv)
if smiles_col not in df.columns:
    raise ValueError(f"SMILES column not found: {smiles_col}")
if smiles_col != "SMILES":
    df = df.rename(columns={smiles_col: "SMILES"})

df = df[df["SMILES"].notna()].copy()
df["SMILES"] = df["SMILES"].astype(str).str.strip()
df = df[df["SMILES"] != ""]
df = df.drop_duplicates(subset=["SMILES"]).reset_index(drop=True)
if df.empty:
    raise ValueError("No compounds available after SMILES normalization.")

df.to_csv(output_csv, index=False)
print(output_csv)
PY
}

prepare_rosettafold_prefold_input_for_protein() {
    local protein=$1
    local rf_dir="${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold/protein_folding"
    local input_dir="${rf_dir}/input"
    local config_dir="${rf_dir}/config"
    local output_dir="${rf_dir}/output"
    local fasta_file="${input_dir}/${protein}.fasta"
    local config_file="${config_dir}/protein_fold.yaml"

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would write RoseTTAFold prefold FASTA/config for ${protein}"
        return 0
    fi

    mkdir -p "$input_dir" "$config_dir" "$output_dir"

    python3 - "$SEQS_CSV" "$protein" "$fasta_file" <<'PY'
import pandas as pd
import sys

seq_csv, protein, fasta_out = sys.argv[1:4]
df = pd.read_csv(seq_csv)
row = df[df["name"] == protein]
if row.empty:
    raise ValueError(f"Protein not found in sequences.csv: {protein}")

sequence = None
for col in ["sequence", "fasta", "seq", "protein_sequence"]:
    if col in row.columns:
        sequence = row[col].iloc[0]
        break
if sequence is None:
    raise ValueError("No supported sequence column found in sequences.csv")

with open(fasta_out, "w") as f:
    f.write(f">{protein}\n{sequence}\n")
PY

    cat > "$config_file" <<EOF
defaults:
  - base

job_name: "${protein}_fold"
output_path: "${output_dir}"

protein_inputs:
  A:
    fasta_file: "${fasta_file}"
EOF
}

stage_1_write_prefold_inputs() {
    log_info "Preparing AF3-first selected_all compound files + RoseTTAFold prefold inputs"
    for protein in $PROTEINS; do
        prepare_all_compounds_for_protein "$protein"
        prepare_rosettafold_prefold_input_for_protein "$protein"
    done

    run_or_dry bash "$SCRIPTS_ROOT/prepare_input/run_write_prefold_af3.sh"
    run_or_dry bash "$SCRIPTS_ROOT/prepare_input/run_write_prefold_boltz2.sh"

    if [ "$DRY_RUN" -eq 1 ]; then return 0; fi

    local ok=1
    for protein in $PROTEINS; do
        local all_csv="${TASK_ROOT}/${protein}/${SELECTED_REL_ALL}"
        local rf_fasta="${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold/protein_folding/input/${protein}.fasta"
        local rf_cfg="${TASK_ROOT}/${protein}/fine_screening/RoseTTAFold/protein_folding/config/protein_fold.yaml"
        if [ ! -f "$all_csv" ]; then
            log_error "Missing selected_all CSV for ${protein}: ${all_csv}"
            ok=0
        fi
        if [ ! -f "$rf_fasta" ]; then
            log_error "Missing RoseTTAFold prefold FASTA for ${protein}: ${rf_fasta}"
            ok=0
        fi
        if [ ! -f "$rf_cfg" ]; then
            log_error "Missing RoseTTAFold prefold config for ${protein}: ${rf_cfg}"
            ok=0
        fi
    done
    [ "$ok" -eq 1 ] || return 1
}

stage_2_prefold_and_vina_receptor() {
    run_or_dry bash "$SCRIPTS_ROOT/structure_prediction/run_prefold_af3.sh"
    run_or_dry bash "$SCRIPTS_ROOT/structure_prediction/run_prefold_boltz2.sh"
    run_or_dry bash "$SCRIPTS_ROOT/structure_prediction/run_rosettafold_prefold.sh"

    for protein in $PROTEINS; do
        local receptor_pdbqt="${TASK_ROOT}/${protein}/fine_screening/Vina/receptor/${protein}.pdbqt"
        run_or_dry bash "$SCRIPTS_ROOT/prepare_input/prepare_vina_receptor.sh" \
            --task-root "$TASK_ROOT" \
            --protein "$protein" \
            --output "$receptor_pdbqt" \
            --conda-env "$VINA_CONDA_ENV"
    done

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would wait for AF3/Boltz2/RoseTTAFold prefold jobs"
        return 0
    fi

    local tmp_jobs
    tmp_jobs=$(mktemp)

    local jf
    jf=$(collect_all_job_files "fine_screening/AF3/prefold")
    [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
    rm -f "$jf"

    jf=$(collect_all_job_files "fine_screening/Boltz2/prefold")
    [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
    rm -f "$jf"

    jf=$(collect_all_job_files "fine_screening/RoseTTAFold")
    [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
    rm -f "$jf"

    wait_for_slurm_jobs "$tmp_jobs" "prefold"
    rm -f "$tmp_jobs"
}

stage_3_run_af3_all_compounds() {
    run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_ALL" bash "$SCRIPTS_ROOT/prepare_input/run_write_af3_input.sh"
    run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_ALL" bash "$SCRIPTS_ROOT/structure_prediction/run_af3_batch.sh"

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would wait for AF3 batch jobs"
        return 0
    fi

    local af3_jobs
    af3_jobs=$(collect_all_job_files "fine_screening/AF3/output")
    wait_for_slurm_jobs "$af3_jobs" "AF3"
    rm -f "$af3_jobs"
}

stage_4_collect_af3_and_select() {
    for protein in $PROTEINS; do
        local base="${TASK_ROOT}/${protein}"
        local af3_dir="${base}/fine_screening/AF3"
        local af3_out="${af3_dir}/output"
        local af3_summary="${af3_dir}/summary.csv"
        local selected_all="${base}/${SELECTED_REL_ALL}"
        local selected_af3="${base}/${SELECTED_REL_AF3}"
        local metrics_csv="${base}/fine_screening/AF3_first/af3_selection_metrics.csv"

        if [ "$DRY_RUN" -eq 1 ]; then
            echo "  [DRY-RUN] Would extract AF3 archives and score/select for ${protein}"
            continue
        fi

        for f in "${af3_out}"/batch_*.tar.gz; do
            [ -f "$f" ] && tar -xzf "$f" -C "${af3_out}/"
        done

        python "$SCRIPTS_ROOT/scoring/af3_scores.py" \
            --af3-results-folder "$af3_out" \
            --output-dir "$af3_dir"

        python "$SCRIPTS_ROOT/scoring/select_af3_top.py" \
            --af3-summary "$af3_summary" \
            --selected-all "$selected_all" \
            --output "$selected_af3" \
            --metrics-output "$metrics_csv" \
            --smiles-col "SMILES" \
            --plddt-threshold 70 \
            --top-fraction 0.5

        if [ ! -f "$selected_af3" ]; then
            log_error "AF3 selected file missing for ${protein}: ${selected_af3}"
            return 1
        fi

        local n_sel
        n_sel=$(tail -n +2 "$selected_af3" | wc -l)
        if [ "$n_sel" -eq 0 ]; then
            log_error "AF3 selected file is empty for ${protein}: ${selected_af3}"
            return 1
        fi

        echo "  ${protein}: AF3 selected ${n_sel} compounds"
    done
}

clear_downstream_outputs() {
    local protein=$1
    local base="${TASK_ROOT}/${protein}"

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would clear downstream outputs for ${protein}"
        return 0
    fi

    rm -rf "${base}/fine_screening/Vina/input" "${base}/fine_screening/Vina/output"
    rm -rf "${base}/fine_screening/Boltz2/output"
    rm -rf "${base}/fine_screening/RoseTTAFold/protein_ligand/output" \
           "${base}/fine_screening/RoseTTAFold/protein_ligand/logs" \
           "${base}/fine_screening/RoseTTAFold/protein_ligand/config"

    rm -f "/home/ubuntu/${protein}/boltz2_tmp/boltz_input.done"
    rm -rf "/home/ubuntu/${protein}/boltz2_tmp/input"
}

stage_5_run_selected_methods() {
    for protein in $PROTEINS; do
        local selected_af3="${TASK_ROOT}/${protein}/${SELECTED_REL_AF3}"
        if [ ! -f "$selected_af3" ] && [ "$DRY_RUN" -eq 0 ]; then
            log_error "Missing AF3-selected compounds for ${protein}: ${selected_af3}"
            return 1
        fi
        clear_downstream_outputs "$protein"
    done

    run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_AF3" bash "$SCRIPTS_ROOT/prepare_input/run_write_boltz2_input.sh"
    run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_AF3" bash "$SCRIPTS_ROOT/docking/run_vina_batch.sh"
    run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_AF3" bash "$SCRIPTS_ROOT/structure_prediction/run_boltz2_batch.sh"
    run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_AF3" bash "$SCRIPTS_ROOT/structure_prediction/run_rosettafold_batch.sh"

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would wait for Vina/Boltz2/RoseTTAFold jobs"
        return 0
    fi

    local tmp_jobs
    tmp_jobs=$(mktemp)

    local jf
    jf=$(collect_all_job_files "fine_screening/Vina/output")
    [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
    rm -f "$jf"

    jf=$(collect_all_job_files "fine_screening/Boltz2/output")
    [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
    rm -f "$jf"

    jf=$(collect_all_job_files "fine_screening/RoseTTAFold")
    [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
    rm -f "$jf"

    wait_for_slurm_jobs "$tmp_jobs" "selected_fine_screening"
    rm -f "$tmp_jobs"
}

stage_6_collect_results() {
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would collect AF3/Vina/Boltz2/RoseTTAFold summaries"
        return 0
    fi

    for protein in $PROTEINS; do
        local base="${TASK_ROOT}/${protein}"
        local selected_af3="${base}/${SELECTED_REL_AF3}"

        # Vina
        python "$SCRIPTS_ROOT/scoring/vina_scores.py" \
            --vina-results-folder "${base}/fine_screening/Vina/output" \
            --output-dir "${base}/fine_screening/Vina" \
            --input-dir "${base}/fine_screening/Vina/input"

        # Boltz2
        local boltz_out="${base}/fine_screening/Boltz2/output"
        for f in "${boltz_out}"/batch_*.tar.gz; do
            [ -f "$f" ] && tar -xzf "$f" -C "${boltz_out}/"
        done
        python "$SCRIPTS_ROOT/scoring/boltz2_scores.py" \
            --boltz-results-folder "$boltz_out" \
            --output-dir "${base}/fine_screening/Boltz2" \
            --smiles "$selected_af3"

        # RoseTTAFold
        python "$SCRIPTS_ROOT/scoring/rosettafold_scores.py" \
            --rfaa-results-folder "${base}/fine_screening/RoseTTAFold/protein_ligand/output" \
            --protein-name "$protein" \
            --output-dir "${base}/fine_screening/RoseTTAFold"
    done
}

show_status() {
    echo "=========================================="
    echo "  AF3-first Pipeline Status"
    echo "  $(date)"
    echo "=========================================="
    echo "Task Root: $TASK_ROOT"
    echo "Proteins:  $PROTEINS"
    echo ""

    for protein in $PROTEINS; do
        local base="${TASK_ROOT}/${protein}"
        echo "=========================================="
        echo "Protein: $protein"
        echo "=========================================="

        local all_csv="${base}/${SELECTED_REL_ALL}"
        local af3_csv="${base}/${SELECTED_REL_AF3}"
        [ -f "$all_csv" ] && echo "  selected_all: present" || echo "  selected_all: missing"
        [ -f "$af3_csv" ] && echo "  selected_af3: present" || echo "  selected_af3: missing"

        [ -f "${base}/fine_screening/AF3/prefold/prefold.done" ] \
            && echo "  AF3 prefold: done" || echo "  AF3 prefold: pending"
        [ -f "/home/ubuntu/${protein}/boltz2_tmp/boltz_results_${protein}/predictions/${protein}/confidence_${protein}_model_0.json" ] \
            && echo "  Boltz2 prefold: done" || echo "  Boltz2 prefold: pending"
        [ -f "${base}/fine_screening/RoseTTAFold/protein_folding/output/protein_fold.done" ] \
            && echo "  RoseTTAFold prefold: done" || echo "  RoseTTAFold prefold: pending"
        [ -f "${base}/fine_screening/Vina/receptor/${protein}.pdbqt" ] \
            && echo "  Vina receptor: ready" || echo "  Vina receptor: missing"

        [ -f "${base}/fine_screening/AF3/summary.csv" ] \
            && echo "  AF3 summary: present" || echo "  AF3 summary: missing"
        [ -f "${base}/fine_screening/Vina/results.csv" ] \
            && echo "  Vina results: present" || echo "  Vina results: missing"
        [ -f "${base}/fine_screening/Boltz2/summary.csv" ] \
            && echo "  Boltz2 summary: present" || echo "  Boltz2 summary: missing"
        [ -f "${base}/fine_screening/RoseTTAFold/summary.csv" ] \
            && echo "  RoseTTAFold summary: present" || echo "  RoseTTAFold summary: missing"

        echo ""
    done
}

show_help() {
    cat <<'HELPEOF'
AF3-first Fine-Screening Pipeline
=================================

Usage: ./scripts/orchestration/run_af3_first_pipeline.sh [OPTIONS] [COMMAND]

Commands:
  all              Run full AF3-first pipeline (default)
  status           Show stage-level status
  help             Show this help message

Options:
  --task-root DIR          Task root directory
  --proteins "P1 P2"       Space-separated protein list
  --smiles-col NAME        SMILES column in Input/compounds_smiles.csv (default: SMILES)
  --start-from STAGE       Resume from stage number 1-6 (default: 1)
  --stop-after STAGE       Stop after stage number 1-6 (default: 6)
  --poll-interval SEC      SLURM poll interval in seconds (default: 300)
  --dry-run                Show what would run without executing

Stages:
  1  write_prefold_inputs      Prepare selected_all.csv + AF3/RoseTTAFold/Boltz2 prefold inputs
  2  prefold_and_receptor      Run AF3/Boltz2/RoseTTAFold prefold + Vina receptor prep
  3  run_af3_all               Run AF3 with all compounds
  4  select_from_af3           Collect AF3 and produce AF3-selected compounds
  5  run_selected_methods      Run Vina/Boltz2/RoseTTAFold on AF3-selected compounds
  6  collect_results           Collect Vina/Boltz2/RoseTTAFold summaries

HELPEOF
}

COMMAND="all"

while [[ $# -gt 0 ]]; do
    case $1 in
        --task-root)
            TASK_ROOT="$2"
            COMPOUNDS_CSV="${TASK_ROOT}/Input/compounds_smiles.csv"
            SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"
            shift 2
            ;;
        --proteins)
            PROTEINS="$2"
            shift 2
            ;;
        --smiles-col)
            SMILES_COL="$2"
            shift 2
            ;;
        --start-from)
            START_FROM="$2"
            shift 2
            ;;
        --stop-after)
            STOP_AFTER="$2"
            shift 2
            ;;
        --poll-interval)
            POLL_INTERVAL="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        all|status|help)
            COMMAND="$1"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run '$SELF_SCRIPT help' for usage information"
            exit 1
            ;;
    esac
done

case $COMMAND in
    help)
        show_help
        exit 0
        ;;
    status)
        show_status
        exit 0
        ;;
esac

if ! preflight_check_scripts; then
    log_error "Preflight script path validation failed"
    exit 1
fi

PIPELINE_LOG="${TASK_ROOT}/af3_first_pipeline.log"
mkdir -p "$TASK_ROOT"

export MASTER_TASK_ROOT="$TASK_ROOT"
export MASTER_PROTEINS="$PROTEINS"

log_info "AF3-first pipeline started"
echo "Task Root:      $TASK_ROOT"
echo "Proteins:       $PROTEINS"
echo "Stages:         ${START_FROM} -> ${STOP_AFTER}"
echo "Poll interval:  ${POLL_INTERVAL}s"
echo "Dry run:        $([ "$DRY_RUN" -eq 1 ] && echo "yes" || echo "no")"
echo ""

run_stage 1 "write_prefold_inputs" stage_1_write_prefold_inputs
run_stage 2 "prefold_and_receptor" stage_2_prefold_and_vina_receptor
run_stage 3 "run_af3_all" stage_3_run_af3_all_compounds
run_stage 4 "select_from_af3" stage_4_collect_af3_and_select
run_stage 5 "run_selected_methods" stage_5_run_selected_methods
run_stage 6 "collect_results" stage_6_collect_results

log_info "AF3-first pipeline finished"
echo "To check results: $SELF_SCRIPT --task-root \"$TASK_ROOT\" --proteins \"$PROTEINS\" status"
echo ""
