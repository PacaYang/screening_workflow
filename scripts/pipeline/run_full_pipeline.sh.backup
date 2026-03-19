#!/bin/bash
#
# Full Automation Pipeline for Screening Workflow
# Orchestrates all 6 stages: make_input → initial_screening → compile_summary
#   → prepare_fine → fine_screening → collect_results
#
# Each stage checks for completion before proceeding, making re-runs safe.
#

set -e

# ============================================================================
# Configuration
# ============================================================================

TASK_ROOT="${MASTER_TASK_ROOT:-/home/yangl_pacagen_com/snake_test}"
PROTEINS="${MASTER_PROTEINS:-JAK1JH1}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Stage control
START_FROM=1
STOP_AFTER=6
POLL_INTERVAL=300
DRY_RUN=0
CURRENT_STAGE=0

# Skip flags — initial screening
SKIP_GRAPHDTA=0
SKIP_HMSA=0
SKIP_COLDDTA=0
SKIP_DRUGLAMP=0
SKIP_CONPLEX=0

# Skip flags — fine screening
SKIP_AF3=0
SKIP_BOLTZ2=0
SKIP_ROSETTAFOLD=0
SKIP_VINA=0
SKIP_DIFFDOCK=0
SKIP_MD_PBSA=0

# Model weights configuration
MODEL_WEIGHTS_DIR="/home/yangl_pacagen_com/Applications/model_weights"

GRAPHDTA_MODEL="${MODEL_WEIGHTS_DIR}/GraphDTA/model_GINConvNet_kiba.pt"
HMSA_MODEL="${MODEL_WEIGHTS_DIR}/HMSA/model.pt"
COLDDTA_CHECKPOINT="${MODEL_WEIGHTS_DIR}/ColdDTA/epoch1297test_loss0.1798.pt"
DRUGLAMP_CHECKPOINT="${MODEL_WEIGHTS_DIR}/DrugLAMP/max_val_ausum= 1.84364.ckpt"
CONPLEX_MODEL="${MODEL_WEIGHTS_DIR}/ConPLex/ConPLex_v1_BindingDB.pt"
AF3_WEIGHT_DIR="${MODEL_WEIGHTS_DIR}/AF3"
AF3_DB_DIR="${MODEL_WEIGHTS_DIR}/af3_db"
RFAA_WEIGHTS="${MODEL_WEIGHTS_DIR}/RoseTTAFold/RFAA_paper_weights.pt"
ROSETTA_DB_UR30="${MODEL_WEIGHTS_DIR}/rosetta_db/UniRef30_2020_06/UniRef30_2020_06"
ROSETTA_DB_BFD="${MODEL_WEIGHTS_DIR}/rosetta_db/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"

# Compile summary configuration
TARGET_N=3000

# ============================================================================
# Logging
# ============================================================================

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

# ============================================================================
# Trap handler
# ============================================================================

cleanup() {
    echo ""
    log_error "Pipeline interrupted at stage ${CURRENT_STAGE}"
    echo "To resume, run:"
    echo "  $0 --start-from ${CURRENT_STAGE} --task-root \"$TASK_ROOT\" --proteins \"$PROTEINS\""
    echo ""
    exit 130
}
trap cleanup SIGINT SIGTERM

# ============================================================================
# Core helpers
# ============================================================================

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

# ============================================================================
# Stage 1: Make Input
# ============================================================================

stage_make_input() {
    run_or_dry bash "$SCRIPT_DIR/../input/run_make_input.sh"

    if [ "$DRY_RUN" -eq 1 ]; then return 0; fi

    # Verify
    local ok=1
    for protein in $PROTEINS; do
        local token="${TASK_ROOT}/${protein}/initial_screening/inputs/finish.token"
        if [ -f "$token" ]; then
            echo "  ${protein}: input ready"
        else
            log_error "${protein}: finish.token not found at ${token}"
            ok=0
        fi
    done
    [ "$ok" -eq 1 ] || { log_error "Stage 1 verification failed"; return 1; }
}

# ============================================================================
# Stage 2: Initial Screening
# ============================================================================

stage_initial_screening() {
    # Submit all enabled methods
    [ "$SKIP_GRAPHDTA" -eq 0 ]  && run_or_dry bash "$SCRIPT_DIR/../initial_screening/run_graphdta.sh"
    [ "$SKIP_HMSA" -eq 0 ]      && run_or_dry bash "$SCRIPT_DIR/../initial_screening/run_hmsa.sh"
    [ "$SKIP_COLDDTA" -eq 0 ]   && run_or_dry bash "$SCRIPT_DIR/../initial_screening/run_colddta.sh"
    [ "$SKIP_DRUGLAMP" -eq 0 ]  && run_or_dry bash "$SCRIPT_DIR/../initial_screening/run_druglamp.sh"
    [ "$SKIP_CONPLEX" -eq 0 ]   && run_or_dry bash "$SCRIPT_DIR/../initial_screening/run_conplex.sh"

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would poll SLURM until all initial screening jobs finish"
        return 0
    fi

    # Collect all job IDs and wait
    local tmp_jobs
    tmp_jobs=$(mktemp)
    local methods="GraphDTA HMSA ColdDTA DrugLAMP ConPLex"
    for method in $methods; do
        local jf
        jf=$(collect_all_job_files "initial_screening/${method}")
        if [ -s "$jf" ]; then cat "$jf" >> "$tmp_jobs"; fi
        rm -f "$jf"
    done
    wait_for_slurm_jobs "$tmp_jobs" "initial_screening"
    rm -f "$tmp_jobs"

    # Verify predictions exist
    local ok=1
    for protein in $PROTEINS; do
        for method in GraphDTA HMSA ColdDTA DrugLAMP ConPLex; do
            case $method in
                GraphDTA)   [ "$SKIP_GRAPHDTA" -eq 1 ] && continue ;;
                HMSA)       [ "$SKIP_HMSA" -eq 1 ] && continue ;;
                ColdDTA)    [ "$SKIP_COLDDTA" -eq 1 ] && continue ;;
                DrugLAMP)   [ "$SKIP_DRUGLAMP" -eq 1 ] && continue ;;
                ConPLex)    [ "$SKIP_CONPLEX" -eq 1 ] && continue ;;
            esac
            local pred="${TASK_ROOT}/${protein}/initial_screening/${method}"
            local found
            found=$(find "$pred" -maxdepth 1 -name "prediction_*.csv" 2>/dev/null | wc -l)
            if [ "$found" -gt 0 ]; then
                echo "  ${protein}/${method}: predictions found"
            else
                log_error "${protein}/${method}: no prediction_*.csv found"
                ok=0
            fi
        done
    done
    [ "$ok" -eq 1 ] || { log_error "Stage 2 verification failed"; return 1; }
}

# ============================================================================
# Stage 3: Compile Summary
# ============================================================================

stage_compile_summary() {
    run_or_dry bash "$SCRIPT_DIR/../scoring/run_compile_summary.sh"

    if [ "$DRY_RUN" -eq 1 ]; then return 0; fi

    local ok=1
    for protein in $PROTEINS; do
        local sel="${TASK_ROOT}/${protein}/initial_screening/selected.csv"
        if [ -f "$sel" ]; then
            local n
            n=$(tail -n +2 "$sel" | wc -l)
            echo "  ${protein}: ${n} compounds selected"
        else
            log_error "${protein}: selected.csv not found"
            ok=0
        fi
    done
    [ "$ok" -eq 1 ] || { log_error "Stage 3 verification failed"; return 1; }
}

# ============================================================================
# Stage 4: Prepare Fine Screening
# ============================================================================

stage_prepare_fine() {
    # 4a — Write prefold inputs
    log_info "Stage 4a: Writing prefold inputs"
    [ "$SKIP_AF3" -eq 0 ]    && run_or_dry bash "$SCRIPT_DIR/../input/run_write_prefold_af3.sh"
    [ "$SKIP_BOLTZ2" -eq 0 ] && run_or_dry bash "$SCRIPT_DIR/../input/run_write_prefold_boltz2.sh"

    # 4b — Submit prefold SLURM jobs
    log_info "Stage 4b: Submitting prefold jobs"
    [ "$SKIP_AF3" -eq 0 ]         && run_or_dry bash "$SCRIPT_DIR/../structure_prediction/prefold/run_prefold_af3.sh"
    [ "$SKIP_BOLTZ2" -eq 0 ]      && run_or_dry bash "$SCRIPT_DIR/../structure_prediction/prefold/run_prefold_boltz2.sh"
    [ "$SKIP_ROSETTAFOLD" -eq 0 ] && run_or_dry bash "$SCRIPT_DIR/../structure_prediction/prefold/run_rosettafold_prefold.sh"

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would poll SLURM until all prefold jobs finish"
        echo "  [DRY-RUN] Would then write fine screening inputs"
        return 0
    fi

    # 4c — Wait for prefold jobs
    log_info "Stage 4c: Waiting for prefold jobs"
    local tmp_jobs
    tmp_jobs=$(mktemp)

    if [ "$SKIP_AF3" -eq 0 ]; then
        local jf
        jf=$(collect_all_job_files "fine_screening/AF3/prefold")
        [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
        rm -f "$jf"
    fi
    if [ "$SKIP_BOLTZ2" -eq 0 ]; then
        local jf
        jf=$(collect_all_job_files "fine_screening/Boltz2/prefold")
        [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
        rm -f "$jf"
    fi
    if [ "$SKIP_ROSETTAFOLD" -eq 0 ]; then
        local jf
        jf=$(collect_all_job_files "fine_screening/RoseTTAFold")
        [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
        rm -f "$jf"
    fi

    wait_for_slurm_jobs "$tmp_jobs" "prefold"
    rm -f "$tmp_jobs"

    # 4d — Write fine screening inputs
    log_info "Stage 4d: Writing fine screening inputs"
    [ "$SKIP_AF3" -eq 0 ]    && bash "$SCRIPT_DIR/../input/run_write_af3_input.sh"
    [ "$SKIP_BOLTZ2" -eq 0 ] && bash "$SCRIPT_DIR/../input/run_write_boltz2_input.sh"
    bash "$SCRIPT_DIR/../input/run_split_csv.sh"
}

# ============================================================================
# Stage 5: Fine Screening
# ============================================================================

stage_fine_screening() {
    # 5a — Submit AF3, Boltz2, Vina, RoseTTAFold
    log_info "Stage 5a: Submitting AF3, Boltz2, Vina, RoseTTAFold jobs"
    [ "$SKIP_AF3" -eq 0 ]         && run_or_dry bash "$SCRIPT_DIR/../structure_prediction/run_af3_batch.sh"
    [ "$SKIP_BOLTZ2" -eq 0 ]      && run_or_dry bash "$SCRIPT_DIR/../structure_prediction/run_boltz2_batch.sh"
    [ "$SKIP_VINA" -eq 0 ]        && run_or_dry bash "$SCRIPT_DIR/../docking/run_vina_batch.sh"
    [ "$SKIP_ROSETTAFOLD" -eq 0 ] && run_or_dry bash "$SCRIPT_DIR/../structure_prediction/run_rosettafold_batch.sh"

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would poll AF3/Boltz2/Vina/RoseTTAFold, then submit DiffDock, poll DiffDock, then submit MD+PBSA, poll MD+PBSA"
        return 0
    fi

    # 5b — Wait for AF3, Boltz2, Vina, RoseTTAFold
    log_info "Stage 5b: Waiting for AF3, Boltz2, Vina, RoseTTAFold jobs"
    local tmp_jobs
    tmp_jobs=$(mktemp)

    local job_paths=""
    [ "$SKIP_AF3" -eq 0 ]         && job_paths="$job_paths fine_screening/AF3/output"
    [ "$SKIP_BOLTZ2" -eq 0 ]      && job_paths="$job_paths fine_screening/Boltz2/output"
    [ "$SKIP_VINA" -eq 0 ]        && job_paths="$job_paths fine_screening/Vina/output"
    [ "$SKIP_ROSETTAFOLD" -eq 0 ] && job_paths="$job_paths fine_screening/RoseTTAFold/protein_ligand"

    for jp in $job_paths; do
        local jf
        jf=$(collect_all_job_files "$jp")
        [ -s "$jf" ] && cat "$jf" >> "$tmp_jobs"
        rm -f "$jf"
    done

    wait_for_slurm_jobs "$tmp_jobs" "AF3/Boltz2/Vina/RoseTTAFold"
    rm -f "$tmp_jobs"

    # 5c — Submit DiffDock
    if [ "$SKIP_DIFFDOCK" -eq 0 ]; then
        log_info "Stage 5c: Submitting DiffDock jobs"
        bash "$SCRIPT_DIR/../docking/run_diffdock_batch.sh"
    fi

    # 5d — Wait for DiffDock
    if [ "$SKIP_DIFFDOCK" -eq 0 ]; then
        log_info "Stage 5d: Waiting for DiffDock jobs"
        local dd_jobs
        dd_jobs=$(collect_all_job_files "fine_screening/PBSA/DiffDock/output")
        wait_for_slurm_jobs "$dd_jobs" "DiffDock"
        rm -f "$dd_jobs"
    fi

    # 5e — Submit MD+PBSA
    if [ "$SKIP_MD_PBSA" -eq 0 ]; then
        log_info "Stage 5e: Submitting MD+PBSA jobs"
        bash "$SCRIPT_DIR/../md_pbsa/run_md_pbsa_batch.sh"
    fi

    # 5f — Wait for MD+PBSA
    if [ "$SKIP_MD_PBSA" -eq 0 ]; then
        log_info "Stage 5f: Waiting for MD+PBSA jobs"
        local pbsa_jobs
        pbsa_jobs=$(collect_all_job_files "fine_screening/PBSA/PBSA")
        wait_for_slurm_jobs "$pbsa_jobs" "MD+PBSA"
        rm -f "$pbsa_jobs"
    fi
}

# ============================================================================
# Stage 6: Collect Results
# ============================================================================

stage_collect_results() {
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "  [DRY-RUN] Would collect scores for AF3, Boltz2, Vina, PBSA per protein"
        return 0
    fi

    for protein in $PROTEINS; do
        log_info "Collecting results for ${protein}"
        local base="${TASK_ROOT}/${protein}"
        local selected="${base}/initial_screening/selected.csv"

        # AF3 scores
        if [ "$SKIP_AF3" -eq 0 ]; then
            local af3_summary="${base}/fine_screening/AF3/summary.csv"
            if [ -f "$af3_summary" ]; then
                echo "  ${protein}/AF3: summary.csv already exists, skipping"
            else
                echo "  ${protein}/AF3: collecting scores"
                local af3_out="${base}/fine_screening/AF3/output"
                python "$SCRIPT_DIR/../scoring/af3_scores.py" \
                    --af3-results-folder "${af3_out}" \
                    --output-dir "${base}/fine_screening/AF3"
            fi
        fi

        # Boltz2 scores
        if [ "$SKIP_BOLTZ2" -eq 0 ]; then
            local boltz2_summary="${base}/fine_screening/Boltz2/summary.csv"
            if [ -f "$boltz2_summary" ]; then
                echo "  ${protein}/Boltz2: summary.csv already exists, skipping"
            else
                echo "  ${protein}/Boltz2: collecting scores"
                local boltz_out="${base}/fine_screening/Boltz2/output"
                python "$SCRIPT_DIR/../scoring/boltz2_scores.py" \
                    --boltz-results-folder "${boltz_out}" \
                    --output-dir "${base}/fine_screening/Boltz2"
            fi
        fi

        # Vina scores
        if [ "$SKIP_VINA" -eq 0 ]; then
            local vina_results="${base}/fine_screening/Vina/results.csv"
            if [ -f "$vina_results" ]; then
                echo "  ${protein}/Vina: results.csv already exists, skipping"
            else
                echo "  ${protein}/Vina: collecting scores"
                python "$SCRIPT_DIR/../scoring/vina_scores.py" \
                    --vina-results-folder "${base}/fine_screening/Vina/output" \
                    --output-dir "${base}/fine_screening/Vina" \
                    --input-dir "${base}/fine_screening/Vina/input"
            fi
        fi

        # PBSA scores
        if [ "$SKIP_MD_PBSA" -eq 0 ]; then
            local pbsa_summary="${base}/fine_screening/PBSA/summary.csv"
            if [ -f "$pbsa_summary" ]; then
                echo "  ${protein}/PBSA: summary.csv already exists, skipping"
            else
                local pbsa_dir="${base}/fine_screening/PBSA/PBSA/PBSA"
                local pbsa_outdir="${base}/fine_screening/PBSA"
                echo "  ${protein}/PBSA: extracting results"
                bash "$SCRIPT_DIR/../md_pbsa/pbsa/pbsa_extract_results.sh" "$pbsa_dir" "$pbsa_outdir"
                echo "  ${protein}/PBSA: mapping SMILES"
                python "$SCRIPT_DIR/../md_pbsa/pbsa/mapping_smiles.py" \
                    --collected "${pbsa_outdir}/tmp.csv" \
                    --smiles_csv "$selected" \
                    --outdir "$pbsa_outdir"
            fi
        fi
    done
}

# ============================================================================
# Status command
# ============================================================================

show_status() {
    echo "=========================================="
    echo "  Full Pipeline Status Report"
    echo "  $(date)"
    echo "=========================================="
    echo ""
    echo "Task Root: $TASK_ROOT"
    echo "Proteins:  $PROTEINS"
    echo ""

    for protein in $PROTEINS; do
        echo "=========================================="
        echo "Protein: $protein"
        echo "=========================================="
        local base="${TASK_ROOT}/${protein}"

        # Stage 1
        echo ""
        echo "  Stage 1 — Make Input:"
        local token="${base}/initial_screening/inputs/finish.token"
        if [ -f "$token" ]; then
            local n_inputs
            n_inputs=$(find "${base}/initial_screening/inputs" -name "*.csv" 2>/dev/null | wc -l)
            echo -e "    ${GREEN}Done${NC} (${n_inputs} input files)"
        else
            echo -e "    ${YELLOW}Not done${NC}"
        fi

        # Stage 2
        echo ""
        echo "  Stage 2 — Initial Screening:"
        for method in GraphDTA HMSA ColdDTA DrugLAMP ConPLex; do
            local pred_count
            pred_count=$(find "${base}/initial_screening/${method}" -maxdepth 1 -name "prediction_*.csv" 2>/dev/null | wc -l)
            if [ "$pred_count" -gt 0 ]; then
                echo -e "    ${method}: ${GREEN}${pred_count} prediction file(s)${NC}"
            else
                echo -e "    ${method}: ${YELLOW}no predictions${NC}"
            fi
        done

        # Stage 3
        echo ""
        echo "  Stage 3 — Compile Summary:"
        local sel="${base}/initial_screening/selected.csv"
        if [ -f "$sel" ]; then
            local n_sel
            n_sel=$(tail -n +2 "$sel" | wc -l)
            echo -e "    ${GREEN}Done${NC} — ${n_sel} compounds selected"
        else
            echo -e "    ${YELLOW}Not done${NC}"
        fi

        # Stage 4
        echo ""
        echo "  Stage 4 — Prepare Fine Screening:"
        [ -f "${base}/fine_screening/AF3/prefold/prefold.done" ] \
            && echo -e "    AF3 prefold: ${GREEN}done${NC}" \
            || echo -e "    AF3 prefold: ${YELLOW}pending${NC}"
        local boltz_conf
        boltz_conf=$(find "${base}/fine_screening/Boltz2/prefold" -name "confidence_*.json" 2>/dev/null | wc -l)
        echo "    Boltz2 prefold: ${boltz_conf} confidence JSON(s)"
        [ -f "${base}/fine_screening/RoseTTAFold/protein_fold.done" ] \
            && echo -e "    RoseTTAFold prefold: ${GREEN}done${NC}" \
            || echo -e "    RoseTTAFold prefold: ${YELLOW}pending${NC}"
        [ -f "${base}/fine_screening/AF3/af3_input.done" ] \
            && echo -e "    AF3 inputs: ${GREEN}written${NC}" \
            || echo -e "    AF3 inputs: ${YELLOW}not written${NC}"
        [ -f "${base}/fine_screening/Boltz2/boltz_input.done" ] \
            && echo -e "    Boltz2 inputs: ${GREEN}written${NC}" \
            || echo -e "    Boltz2 inputs: ${YELLOW}not written${NC}"
        local vina_chunks
        vina_chunks=$(find "${base}/fine_screening/Vina/input" -name "*.csv" 2>/dev/null | wc -l)
        echo "    Vina/DiffDock input chunks: ${vina_chunks}"

        # Stage 5
        echo ""
        echo "  Stage 5 — Fine Screening:"
        if [ -f "$sel" ]; then
            local n_compounds
            n_compounds=$(tail -n +2 "$sel" | wc -l)

            local af3_done
            af3_done=$(find "${base}/fine_screening/AF3/output" -name "*_summary_confidences.json" 2>/dev/null | wc -l)
            echo "    AF3: ${af3_done}/${n_compounds} compounds"

            local boltz_batches
            boltz_batches=$(find "${base}/fine_screening/Boltz2/output/token" -name "batch_*.done" 2>/dev/null | wc -l)
            echo "    Boltz2: ${boltz_batches}/40 batches"

            local vina_done
            vina_done=$(find "${base}/fine_screening/Vina/output" -name "*.done" -maxdepth 1 2>/dev/null | wc -l)
            echo "    Vina: ${vina_done}/${vina_chunks} parts"

            local dd_done
            dd_done=$(find "${base}/fine_screening/PBSA/DiffDock/output" -name "*.done" -maxdepth 1 2>/dev/null | wc -l)
            echo "    DiffDock: ${dd_done}/${vina_chunks} parts"

            local md_done
            md_done=$(find "${base}/fine_screening/PBSA/PBSA/MD" -name "token.done" 2>/dev/null | wc -l)
            local pbsa_done
            pbsa_done=$(find "${base}/fine_screening/PBSA/PBSA/PBSA" -name "token.done" 2>/dev/null | wc -l)
            echo "    MD: ${md_done}/${n_compounds}  PBSA: ${pbsa_done}/${n_compounds}"
        else
            echo "    (selected.csv not found — cannot determine progress)"
        fi

        # Stage 6
        echo ""
        echo "  Stage 6 — Collected Results:"
        [ -f "${base}/fine_screening/AF3/summary.csv" ] \
            && echo -e "    AF3: ${GREEN}summary.csv exists${NC}" \
            || echo -e "    AF3: ${YELLOW}not collected${NC}"
        [ -f "${base}/fine_screening/Boltz2/summary.csv" ] \
            && echo -e "    Boltz2: ${GREEN}summary.csv exists${NC}" \
            || echo -e "    Boltz2: ${YELLOW}not collected${NC}"
        [ -f "${base}/fine_screening/Vina/results.csv" ] \
            && echo -e "    Vina: ${GREEN}results.csv exists${NC}" \
            || echo -e "    Vina: ${YELLOW}not collected${NC}"
        [ -f "${base}/fine_screening/PBSA/summary.csv" ] \
            && echo -e "    PBSA: ${GREEN}summary.csv exists${NC}" \
            || echo -e "    PBSA: ${YELLOW}not collected${NC}"

        echo ""
    done
}

# ============================================================================
# Help
# ============================================================================

show_help() {
    cat <<'HELPEOF'
Full Automation Pipeline for Screening Workflow
================================================

Usage: run_full_pipeline.sh [OPTIONS] [COMMAND]

Commands:
  all              Run full pipeline (default)
  status           Show progress across all stages
  help             Show this help message

Options:
  --task-root DIR          Task root directory
  --proteins "P1 P2"       Space-separated protein list
  --start-from STAGE       Resume from stage number 1-6 (default: 1)
  --stop-after STAGE       Stop after stage number 1-6 (default: 6)
  --poll-interval SEC      SLURM poll interval in seconds (default: 300)
  --skip-graphdta          Skip GraphDTA in initial screening
  --skip-hmsa              Skip HMSA
  --skip-colddta           Skip ColdDTA
  --skip-druglamp          Skip DrugLAMP
  --skip-conplex           Skip ConPLex
  --skip-af3               Skip AF3 (prefold + fine)
  --skip-boltz2            Skip Boltz2 (prefold + fine)
  --skip-rosettafold       Skip RoseTTAFold (prefold + fine)
  --skip-vina              Skip Vina
  --skip-diffdock          Skip DiffDock
  --skip-md-pbsa           Skip MD+PBSA
  --target-n N             Number of compounds for fine screening (default: 3000)
  --dry-run                Show what would run without executing
  --model-weights-dir DIR  Base directory for all model weights

Stages:
  1  make_input          Prepare input files (head node)
  2  initial_screening   Run 5 screening methods via SLURM, wait for completion
  3  compile_summary     Compile results and select compounds (head node)
  4  prepare_fine        Write prefold inputs, run prefolds via SLURM, write fine inputs
  5  fine_screening      Submit fine screening SLURM jobs, wait for completion
  6  collect_results     Collect and summarize scores (head node)

Examples:
  # Run the full pipeline end-to-end
  run_full_pipeline.sh --task-root /data/screen --proteins "JAK1JH1 EGFR"

  # Resume from stage 4 after fixing a prefold issue
  run_full_pipeline.sh --start-from 4 --task-root /data/screen --proteins "JAK1JH1"

  # Run only initial screening (stages 1-3)
  run_full_pipeline.sh --stop-after 3

  # Dry run to see what would execute
  run_full_pipeline.sh --dry-run all

  # Check progress
  run_full_pipeline.sh status

HELPEOF
}

# ============================================================================
# CLI Argument Parsing
# ============================================================================

COMMAND="all"

while [[ $# -gt 0 ]]; do
    case $1 in
        --task-root)
            TASK_ROOT="$2"
            shift 2
            ;;
        --proteins)
            PROTEINS="$2"
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
        --skip-graphdta)
            SKIP_GRAPHDTA=1
            shift
            ;;
        --skip-hmsa)
            SKIP_HMSA=1
            shift
            ;;
        --skip-colddta)
            SKIP_COLDDTA=1
            shift
            ;;
        --skip-druglamp)
            SKIP_DRUGLAMP=1
            shift
            ;;
        --skip-conplex)
            SKIP_CONPLEX=1
            shift
            ;;
        --skip-af3)
            SKIP_AF3=1
            shift
            ;;
        --skip-boltz2)
            SKIP_BOLTZ2=1
            shift
            ;;
        --skip-rosettafold)
            SKIP_ROSETTAFOLD=1
            shift
            ;;
        --skip-vina)
            SKIP_VINA=1
            shift
            ;;
        --skip-diffdock)
            SKIP_DIFFDOCK=1
            shift
            ;;
        --skip-md-pbsa)
            SKIP_MD_PBSA=1
            shift
            ;;
        --target-n)
            TARGET_N="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --model-weights-dir)
            MODEL_WEIGHTS_DIR="$2"
            GRAPHDTA_MODEL="${MODEL_WEIGHTS_DIR}/GraphDTA/model_GINConvNet_kiba.pt"
            HMSA_MODEL="${MODEL_WEIGHTS_DIR}/HMSA/model.pt"
            COLDDTA_CHECKPOINT="${MODEL_WEIGHTS_DIR}/ColdDTA/epoch1297test_loss0.1798.pt"
            DRUGLAMP_CHECKPOINT="${MODEL_WEIGHTS_DIR}/DrugLAMP/max_val_ausum= 1.84364.ckpt"
            CONPLEX_MODEL="${MODEL_WEIGHTS_DIR}/ConPLex/ConPLex_v1_BindingDB.pt"
            AF3_WEIGHT_DIR="${MODEL_WEIGHTS_DIR}/AF3"
            AF3_DB_DIR="${MODEL_WEIGHTS_DIR}/af3_db"
            RFAA_WEIGHTS="${MODEL_WEIGHTS_DIR}/RoseTTAFold/RFAA_paper_weights.pt"
            ROSETTA_DB_UR30="${MODEL_WEIGHTS_DIR}/rosetta_db/UniRef30_2020_06/UniRef30_2020_06"
            ROSETTA_DB_BFD="${MODEL_WEIGHTS_DIR}/rosetta_db/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
            shift 2
            ;;
        all|status|help)
            COMMAND="$1"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run '$0 help' for usage information"
            exit 1
            ;;
    esac
done

# ============================================================================
# Main
# ============================================================================

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

# Setup pipeline log
PIPELINE_LOG="${TASK_ROOT}/pipeline.log"
mkdir -p "$TASK_ROOT"

# Export config for sub-scripts
export MASTER_TASK_ROOT="$TASK_ROOT"
export MASTER_PROTEINS="$PROTEINS"
export MASTER_GRAPHDTA_MODEL="$GRAPHDTA_MODEL"
export MASTER_HMSA_MODEL="$HMSA_MODEL"
export MASTER_COLDDTA_CHECKPOINT="$COLDDTA_CHECKPOINT"
export MASTER_DRUGLAMP_CHECKPOINT="$DRUGLAMP_CHECKPOINT"
export MASTER_CONPLEX_MODEL="$CONPLEX_MODEL"
export MASTER_AF3_WEIGHT_DIR="$AF3_WEIGHT_DIR"
export MASTER_AF3_DB_DIR="$AF3_DB_DIR"
export MASTER_RFAA_WEIGHTS="$RFAA_WEIGHTS"
export MASTER_ROSETTA_DB_UR30="$ROSETTA_DB_UR30"
export MASTER_ROSETTA_DB_BFD="$ROSETTA_DB_BFD"
export MASTER_TARGET_N="$TARGET_N"

# Banner
echo "=========================================="
echo "  Full Screening Pipeline"
echo "=========================================="
echo ""
echo "Task Root:      $TASK_ROOT"
echo "Proteins:       $PROTEINS"
echo "Weights dir:    $MODEL_WEIGHTS_DIR"
echo "Target N:       $TARGET_N"
echo "Stages:         ${START_FROM} -> ${STOP_AFTER}"
echo "Poll interval:  ${POLL_INTERVAL}s"
echo "Dry run:        $([ "$DRY_RUN" -eq 1 ] && echo "yes" || echo "no")"
echo ""

log_info "Pipeline started"

run_stage 1 "make_input"        stage_make_input
run_stage 2 "initial_screening" stage_initial_screening
run_stage 3 "compile_summary"   stage_compile_summary
run_stage 4 "prepare_fine"      stage_prepare_fine
run_stage 5 "fine_screening"    stage_fine_screening
run_stage 6 "collect_results"   stage_collect_results

log_info "Pipeline finished"
echo "To check results: $0 --task-root \"$TASK_ROOT\" --proteins \"$PROTEINS\" status"
echo ""