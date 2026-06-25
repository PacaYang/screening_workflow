#!/bin/bash
#
# TL1A Fine-Screening Orchestrator
# Runs AF3, Boltz2, RFAA, and Vina on the full compound set.
# Assumes prefold (AF3, Boltz2, RFAA) is already complete.
#

set -e

TASK_ROOT="${MASTER_TASK_ROOT:-/shared/B3/TL1A_Jun10}"
PROTEINS="${MASTER_PROTEINS:-TL1A}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SELF_SCRIPT="${SCRIPTS_ROOT}/orchestration/run_fineScreening_TL1A.sh"

START_FROM=1
STOP_AFTER=4
POLL_INTERVAL=300
DRY_RUN=0
CURRENT_STAGE=0

COMPOUNDS_CSV="${TASK_ROOT}/Input/compounds_smiles.csv"
SEQS_CSV="${TASK_ROOT}/Input/sequences.csv"
SMILES_COL="SMILES"
VINA_CONDA_ENV="${VINA_CONDA_ENV:-vina_new}"

SKIP_AF3=0
SKIP_BOLTZ2=0
SKIP_RFAA=0
SKIP_VINA=0

SELECTED_REL_ALL="fine_screening/AF3_first/selected_all.csv"

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
    echo "$msg" >&2
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
        "$SCRIPTS_ROOT/prepare_input/prepare_vina_receptor.sh"
        "$SCRIPTS_ROOT/prepare_input/run_write_af3_input.sh"
        "$SCRIPTS_ROOT/prepare_input/run_write_boltz2_input.sh"
        "$SCRIPTS_ROOT/structure_prediction/run_af3_batch.sh"
        "$SCRIPTS_ROOT/structure_prediction/run_boltz2_batch.sh"
        "$SCRIPTS_ROOT/structure_prediction/run_rosettafold_batch.sh"
        "$SCRIPTS_ROOT/docking/run_vina_batch.sh"
        "$SCRIPTS_ROOT/scoring/af3_scores.py"
        "$SCRIPTS_ROOT/scoring/boltz2_scores.py"
        "$SCRIPTS_ROOT/scoring/rosettafold_scores.py"
        "$SCRIPTS_ROOT/scoring/vina_scores.py"
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
        echo "  [DRY-RUN] Would write all-compounds CSV: ${out_csv}"
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

stage_1_write_selected_all() {
    log_info "Writing selected_all.csv from compounds_smiles.csv"
    for protein in $PROTEINS; do
        prepare_all_compounds_for_protein "$protein"
    done

    if [ "$DRY_RUN" -eq 1 ]; then return 0; fi

    local ok=1
    for protein in $PROTEINS; do
        local all_csv="${TASK_ROOT}/${protein}/${SELECTED_REL_ALL}"
        if [ ! -f "$all_csv" ]; then
            log_error "Missing selected_all CSV for ${protein}: ${all_csv}"
            ok=0
        fi
    done
    [ "$ok" -eq 1 ] || return 1
}

stage_2_vina_receptor() {
    for protein in $PROTEINS; do
        local receptor_pdbqt="${TASK_ROOT}/${protein}/fine_screening/Vina/receptor/${protein}.pdbqt"
        run_or_dry bash "$SCRIPTS_ROOT/prepare_input/prepare_vina_receptor.sh" \
            --task-root "$TASK_ROOT" \
            --protein "$protein" \
            --output "$receptor_pdbqt" \
            --conda-env "$VINA_CONDA_ENV"
    done
}

stage_3_run_methods() {
    for protein in $PROTEINS; do
        local all_csv="${TASK_ROOT}/${protein}/${SELECTED_REL_ALL}"
        if [ ! -f "$all_csv" ] && [ "$DRY_RUN" -eq 0 ]; then
            log_error "Missing selected_all CSV for ${protein}: ${all_csv}"
            return 1
        fi
    done

    if [ "$SKIP_AF3" -eq 0 ]; then
        run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_ALL" bash "$SCRIPTS_ROOT/prepare_input/run_write_af3_input.sh"
        run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_ALL" bash "$SCRIPTS_ROOT/structure_prediction/run_af3_batch.sh"
    fi

    if [ "$SKIP_BOLTZ2" -eq 0 ]; then
        run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_ALL" bash "$SCRIPTS_ROOT/prepare_input/run_write_boltz2_input.sh"
        run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_ALL" bash "$SCRIPTS_ROOT/structure_prediction/run_boltz2_batch.sh"
    fi

    if [ "$SKIP_RFAA" -eq 0 ]; then
        run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_ALL" bash "$SCRIPTS_ROOT/structure_prediction/run_rosettafold_batch.sh"
    fi

    if [ "$SKIP_VINA" -eq 0 ]; then
        run_or_dry env MASTER_SELECTED_REL_PATH="$SELECTED_REL_ALL" bash "$SCRIPTS_ROOT/docking/run_vina_batch.sh"
    fi

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would wait for AF3/Boltz2/RFAA/Vina jobs"
        return 0
    fi

    local tmp_jobs
    tmp_jobs=$(mktemp)

    local jf
    if [ "$SKIP_AF3" -eq 0 ]; then
        jf=$(collect_all_job_files "fine_screening/AF3/output")
        [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
        rm -f "$jf"
    fi
    if [ "$SKIP_BOLTZ2" -eq 0 ]; then
        jf=$(collect_all_job_files "fine_screening/Boltz2/output")
        [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
        rm -f "$jf"
    fi
    if [ "$SKIP_RFAA" -eq 0 ]; then
        jf=$(collect_all_job_files "fine_screening/RoseTTAFold")
        [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
        rm -f "$jf"
    fi
    if [ "$SKIP_VINA" -eq 0 ]; then
        jf=$(collect_all_job_files "fine_screening/Vina/output")
        [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
        rm -f "$jf"
    fi

    wait_for_slurm_jobs "$tmp_jobs" "fine_screening_4methods"
    rm -f "$tmp_jobs"
}

stage_4_collect_results() {
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would collect AF3/Boltz2/RFAA/Vina summaries"
        return 0
    fi

    for protein in $PROTEINS; do
        local base="${TASK_ROOT}/${protein}"
        local selected_all="${base}/${SELECTED_REL_ALL}"

        if [ "$SKIP_AF3" -eq 0 ]; then
            local af3_dir="${base}/fine_screening/AF3"
            local af3_out="${af3_dir}/output"
            for f in "${af3_out}"/batch_*.tar.gz; do
                [ -f "$f" ] && tar -xzf "$f" -C "${af3_out}/"
            done
            python "$SCRIPTS_ROOT/scoring/af3_scores.py" \
                --af3-results-folder "$af3_out" \
                --output-dir "$af3_dir"
        fi

        if [ "$SKIP_BOLTZ2" -eq 0 ]; then
            local boltz_out="${base}/fine_screening/Boltz2/output"
            for f in "${boltz_out}"/batch_*.tar.gz; do
                [ -f "$f" ] && tar -xzf "$f" -C "${boltz_out}/"
            done
            python "$SCRIPTS_ROOT/scoring/boltz2_scores.py" \
                --boltz-results-folder "$boltz_out" \
                --output-dir "${base}/fine_screening/Boltz2" \
                --smiles "$selected_all"
        fi

        if [ "$SKIP_RFAA" -eq 0 ]; then
            python "$SCRIPTS_ROOT/scoring/rosettafold_scores.py" \
                --rfaa-results-folder "${base}/fine_screening/RoseTTAFold/protein_ligand/output" \
                --protein-name "$protein" \
                --output-dir "${base}/fine_screening/RoseTTAFold"
        fi

        if [ "$SKIP_VINA" -eq 0 ]; then
            python "$SCRIPTS_ROOT/scoring/vina_scores.py" \
                --vina-results-folder "${base}/fine_screening/Vina/output" \
                --output-dir "${base}/fine_screening/Vina" \
                --input-dir "${base}/fine_screening/Vina/input"
        fi
    done
}

show_status() {
    echo "=========================================="
    echo "  TL1A Fine-Screening Status"
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
        [ -f "$all_csv" ] && echo "  selected_all: present ($(tail -n +2 "$all_csv" | wc -l) compounds)" || echo "  selected_all: missing"

        [ -f "${base}/fine_screening/AF3/prefold/prefold.done" ] \
            && echo "  AF3 prefold: done" || echo "  AF3 prefold: pending"
        [ -f "/home/ubuntu/${protein}/boltz2_tmp/boltz_results_${protein}/predictions/${protein}/confidence_${protein}_model_0.json" ] \
            && echo "  Boltz2 prefold: done" || echo "  Boltz2 prefold: pending"
        [ -f "${base}/fine_screening/RoseTTAFold/protein_folding/output/protein_fold.done" ] \
            && echo "  RFAA prefold: done" || echo "  RFAA prefold: pending"
        [ -f "${base}/fine_screening/Vina/receptor/${protein}.pdbqt" ] \
            && echo "  Vina receptor: ready" || echo "  Vina receptor: missing"

        [ -f "${base}/fine_screening/AF3/summary.csv" ] \
            && echo "  AF3 summary: present" || echo "  AF3 summary: missing"
        [ -f "${base}/fine_screening/Boltz2/summary.csv" ] \
            && echo "  Boltz2 summary: present" || echo "  Boltz2 summary: missing"
        [ -f "${base}/fine_screening/RoseTTAFold/summary.csv" ] \
            && echo "  RFAA summary: present" || echo "  RFAA summary: missing"
        [ -f "${base}/fine_screening/Vina/results.csv" ] \
            && echo "  Vina results: present" || echo "  Vina results: missing"

        echo ""
    done
}

show_help() {
    cat <<'HELPEOF'
TL1A Fine-Screening Orchestrator (4 methods: AF3, Boltz2, RFAA, Vina)
=====================================================================

Usage: ./scripts/orchestration/run_fineScreening_TL1A.sh [OPTIONS] [COMMAND]

Commands:
  all              Run full pipeline (default)
  status           Show stage-level status
  help             Show this help message

Options:
  --task-root DIR          Task root directory (default: /shared/B3/TL1A_Jun10)
  --proteins "P1 P2"       Space-separated protein list (default: TL1A)
  --smiles-col NAME        SMILES column in compounds_smiles.csv (default: SMILES)
  --start-from STAGE       Resume from stage 1-4 (default: 1)
  --stop-after STAGE       Stop after stage 1-4 (default: 4)
  --poll-interval SEC      SLURM poll interval in seconds (default: 300)
  --skip-af3               Skip AF3 in stages 3 and 4
  --skip-boltz2            Skip Boltz2 in stages 3 and 4
  --skip-rfaa              Skip RFAA in stages 3 and 4
  --skip-vina              Skip Vina in stages 3 and 4
  --dry-run                Show what would run without executing

Stages:
  1  write_selected_all    Write selected_all.csv from Input/compounds_smiles.csv
  2  vina_receptor         Prepare Vina receptor PDBQT
  3  run_methods           Submit AF3, Boltz2, RFAA, Vina batches and wait
  4  collect_results       Collect AF3/Boltz2/RFAA/Vina summaries

Assumes prefold (AF3, Boltz2, RFAA) is already complete.

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
        --skip-af3)
            SKIP_AF3=1
            shift
            ;;
        --skip-boltz2)
            SKIP_BOLTZ2=1
            shift
            ;;
        --skip-rfaa)
            SKIP_RFAA=1
            shift
            ;;
        --skip-vina)
            SKIP_VINA=1
            shift
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

PIPELINE_LOG="${TASK_ROOT}/fineScreening_TL1A_pipeline.log"
mkdir -p "$TASK_ROOT"

export MASTER_TASK_ROOT="$TASK_ROOT"
export MASTER_PROTEINS="$PROTEINS"

log_info "TL1A fine-screening pipeline started"
echo "Task Root:      $TASK_ROOT"
echo "Proteins:       $PROTEINS"
echo "Stages:         ${START_FROM} -> ${STOP_AFTER}"
echo "Poll interval:  ${POLL_INTERVAL}s"
echo "Methods:        AF3=$([ $SKIP_AF3 -eq 0 ] && echo on || echo off) Boltz2=$([ $SKIP_BOLTZ2 -eq 0 ] && echo on || echo off) RFAA=$([ $SKIP_RFAA -eq 0 ] && echo on || echo off) Vina=$([ $SKIP_VINA -eq 0 ] && echo on || echo off)"
echo "Dry run:        $([ "$DRY_RUN" -eq 1 ] && echo "yes" || echo "no")"
echo ""

run_stage 1 "write_selected_all" stage_1_write_selected_all
run_stage 2 "vina_receptor"      stage_2_vina_receptor
run_stage 3 "run_methods"        stage_3_run_methods
run_stage 4 "collect_results"    stage_4_collect_results

log_info "TL1A fine-screening pipeline finished"
echo "To check results: $SELF_SCRIPT --task-root \"$TASK_ROOT\" --proteins \"$PROTEINS\" status"
echo ""
