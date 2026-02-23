#!/bin/bash
#
# Master Automation Script for Screening Workflow
# This script orchestrates Boltz2, Vina, DiffDock, and MD+PBSA workflows
#

set -e

# ============================================================================
# Configuration - Update these paths as needed
# ============================================================================

# Task root directory
TASK_ROOT="/home/yangl_pacagen_com/snake_test"

# Protein list (space-separated)
PROTEINS="JAK1JH1"

# Script directory (where the automation scripts are located)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Individual workflow scripts
BOLTZ2_SCRIPT="${SCRIPT_DIR}/run_boltz2_batch.sh"
VINA_SCRIPT="${SCRIPT_DIR}/run_vina_batch.sh"
DIFFDOCK_SCRIPT="${SCRIPT_DIR}/run_diffdock_batch.sh"
MD_PBSA_SCRIPT="${SCRIPT_DIR}/run_md_pbsa_batch.sh"

# Workflow control flags (set to 1 to enable, 0 to disable)
RUN_BOLTZ2=1
RUN_VINA=1
RUN_DIFFDOCK=1
RUN_MD_PBSA=1

# Dependency control - wait for DiffDock before MD+PBSA
WAIT_FOR_DIFFDOCK=1  # Set to 1 to wait for DiffDock completion before running MD+PBSA
WAIT_POLL_INTERVAL=300  # Check every 5 minutes (300 seconds)

# ============================================================================
# Functions
# ============================================================================

log_info() {
    echo "=========================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
    echo "=========================================="
    echo ""
}

log_error() {
    echo "=========================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*"
    echo "=========================================="
    echo ""
}

check_script_exists() {
    local script=$1
    if [ ! -f "$script" ]; then
        log_error "Script not found: $script"
        return 1
    fi
    if [ ! -x "$script" ]; then
        log_error "Script not executable: $script"
        log_info "Run: chmod +x $script"
        return 1
    fi
    return 0
}

# Export configuration for sub-scripts
export_config() {
    export MASTER_TASK_ROOT="$TASK_ROOT"
    export MASTER_PROTEINS="$PROTEINS"
    log_info "Configuration exported to environment"
    echo "TASK_ROOT: $TASK_ROOT"
    echo "PROTEINS: $PROTEINS"
    echo ""
}

# ============================================================================
# Workflow Functions
# ============================================================================

run_boltz2() {
    log_info "Starting Boltz2 workflow"

    if ! check_script_exists "$BOLTZ2_SCRIPT"; then
        log_error "Cannot run Boltz2 workflow"
        return 1
    fi

    echo "Running: $BOLTZ2_SCRIPT"
    echo ""

    bash "$BOLTZ2_SCRIPT"

    log_info "Boltz2 workflow submission completed"
    echo "Monitor with: squeue -u \$USER | grep boltz2"
    echo ""
}

run_vina() {
    log_info "Starting Vina workflow"

    if ! check_script_exists "$VINA_SCRIPT"; then
        log_error "Cannot run Vina workflow"
        return 1
    fi

    echo "Running: $VINA_SCRIPT"
    echo ""

    bash "$VINA_SCRIPT"

    log_info "Vina workflow submission completed"
    echo "Monitor with: squeue -u \$USER | grep vina"
    echo ""
}

run_diffdock() {
    log_info "Starting DiffDock workflow"

    if ! check_script_exists "$DIFFDOCK_SCRIPT"; then
        log_error "Cannot run DiffDock workflow"
        return 1
    fi

    echo "Running: $DIFFDOCK_SCRIPT"
    echo ""

    bash "$DIFFDOCK_SCRIPT"

    log_info "DiffDock workflow submission completed"
    echo "Monitor with: squeue -u \$USER | grep diffdock"
    echo ""
}

run_md_pbsa() {
    log_info "Starting MD + PBSA workflow"

    if ! check_script_exists "$MD_PBSA_SCRIPT"; then
        log_error "Cannot run MD+PBSA workflow"
        return 1
    fi

    echo "Running: $MD_PBSA_SCRIPT"
    echo ""

    bash "$MD_PBSA_SCRIPT"

    log_info "MD + PBSA workflow submission completed"
    echo "Monitor with: squeue -u \$USER | grep md_pbsa"
    echo ""
}

# Check if DiffDock is complete for all proteins
check_diffdock_complete() {
    local all_complete=1

    for protein in $PROTEINS; do
        local diffdock_dir="${TASK_ROOT}/${protein}/fine_screening/PBSA/DiffDock/output"
        local input_dir="${TASK_ROOT}/${protein}/fine_screening/Vina/input"

        if [ ! -d "$input_dir" ]; then
            echo "Warning: Input directory not found for ${protein}: $input_dir"
            all_complete=0
            continue
        fi

        # Count expected parts
        local expected_parts=$(ls "$input_dir"/*.csv 2>/dev/null | wc -l)

        if [ $expected_parts -eq 0 ]; then
            echo "Warning: No input parts found for ${protein}"
            all_complete=0
            continue
        fi

        # Count completed parts
        local completed_parts=0
        if [ -d "$diffdock_dir" ]; then
            completed_parts=$(find "$diffdock_dir" -maxdepth 1 -name "*.done" 2>/dev/null | wc -l)
        fi

        echo "${protein}: DiffDock ${completed_parts}/${expected_parts} parts completed"

        if [ $completed_parts -lt $expected_parts ]; then
            all_complete=0
        fi
    done

    return $all_complete
}

# Wait for DiffDock to complete
wait_for_diffdock() {
    log_info "Checking DiffDock completion status"

    if check_diffdock_complete; then
        log_info "DiffDock is already complete for all proteins"
        return 0
    fi

    log_info "DiffDock is not complete. Waiting for completion..."
    echo "Check interval: ${WAIT_POLL_INTERVAL} seconds"
    echo "Press Ctrl+C to stop waiting and skip MD+PBSA"
    echo ""

    local wait_count=0
    while true; do
        sleep $WAIT_POLL_INTERVAL
        wait_count=$((wait_count + 1))

        echo "=========================================="
        echo "Wait check #${wait_count} at $(date)"
        echo "=========================================="

        if check_diffdock_complete; then
            echo ""
            log_info "DiffDock completed! Proceeding to MD+PBSA"
            return 0
        fi

        echo ""
        echo "Still waiting... (checked ${wait_count} times)"
        echo "Next check in ${WAIT_POLL_INTERVAL} seconds"
        echo ""
    done
}

# ============================================================================
# Status Check Functions
# ============================================================================

check_prerequisites() {
    log_info "Checking prerequisites"

    local all_ok=1

    # Check TASK_ROOT exists
    if [ ! -d "$TASK_ROOT" ]; then
        log_error "TASK_ROOT directory not found: $TASK_ROOT"
        all_ok=0
    else
        echo "✓ TASK_ROOT exists: $TASK_ROOT"
    fi

    # Check proteins
    for protein in $PROTEINS; do
        local protein_dir="${TASK_ROOT}/${protein}"
        if [ ! -d "$protein_dir" ]; then
            echo "⚠ Warning: Protein directory not found: $protein_dir"
        else
            echo "✓ Protein directory exists: $protein"
        fi
    done

    echo ""

    # Check for initial screening results
    for protein in $PROTEINS; do
        local selected_csv="${TASK_ROOT}/${protein}/initial_screening/selected.csv"
        if [ ! -f "$selected_csv" ]; then
            echo "⚠ Warning: Initial screening results not found for ${protein}"
            echo "  Expected: $selected_csv"
        else
            local n_compounds=$(tail -n +2 "$selected_csv" | wc -l)
            echo "✓ ${protein}: Found ${n_compounds} selected compounds"
        fi
    done

    echo ""

    if [ $all_ok -eq 0 ]; then
        log_error "Prerequisites check failed"
        return 1
    fi

    log_info "Prerequisites check passed"
    return 0
}

show_status() {
    log_info "Workflow Status"

    for protein in $PROTEINS; do
        echo "Protein: $protein"
        echo "----------------------------------------"

        local protein_dir="${TASK_ROOT}/${protein}"

        if [ ! -d "$protein_dir" ]; then
            echo "  Protein directory not found"
            echo ""
            continue
        fi

        # Boltz2 status
        local boltz2_dir="${protein_dir}/fine_screening/Boltz2/output"
        if [ -d "$boltz2_dir" ]; then
            local boltz2_done=$(find "$boltz2_dir/token" -name "batch_*.done" 2>/dev/null | wc -l)
            echo "  Boltz2: ${boltz2_done}/40 batches completed"
        else
            echo "  Boltz2: Not started"
        fi

        # Vina status
        local vina_dir="${protein_dir}/fine_screening/Vina/output"
        if [ -d "$vina_dir" ]; then
            local vina_done=$(find "$vina_dir" -maxdepth 1 -name "*.done" 2>/dev/null | wc -l)
            local vina_total=$(ls "${protein_dir}/fine_screening/Vina/input"/*.csv 2>/dev/null | wc -l)
            echo "  Vina: ${vina_done}/${vina_total} parts completed"
        else
            echo "  Vina: Not started"
        fi

        # DiffDock status
        local diffdock_dir="${protein_dir}/fine_screening/PBSA/DiffDock/output"
        if [ -d "$diffdock_dir" ]; then
            local diffdock_done=$(find "$diffdock_dir" -maxdepth 1 -name "*.done" 2>/dev/null | wc -l)
            local diffdock_total=$(ls "${protein_dir}/fine_screening/Vina/input"/*.csv 2>/dev/null | wc -l)
            echo "  DiffDock: ${diffdock_done}/${diffdock_total} parts completed"
        else
            echo "  DiffDock: Not started"
        fi

        # MD+PBSA status
        local pbsa_dir="${protein_dir}/fine_screening/PBSA/PBSA"
        if [ -d "${pbsa_dir}/MD" ]; then
            local md_done=$(find "${pbsa_dir}/MD" -mindepth 2 -maxdepth 2 -name "token.done" 2>/dev/null | wc -l)
            local pbsa_done=$(find "${pbsa_dir}/PBSA" -mindepth 2 -maxdepth 2 -name "token.done" 2>/dev/null | wc -l)
            local selected_csv="${protein_dir}/initial_screening/selected.csv"
            local total=0
            if [ -f "$selected_csv" ]; then
                total=$(tail -n +2 "$selected_csv" | wc -l)
            fi
            echo "  MD: ${md_done}/${total} completed"
            echo "  PBSA: ${pbsa_done}/${total} completed"
        else
            echo "  MD+PBSA: Not started"
        fi

        echo ""
    done

    # Show SLURM queue
    echo "SLURM Queue Status:"
    echo "----------------------------------------"
    echo "Running jobs:"
    squeue -u $USER -t R -o "  %.50j %.8T" 2>/dev/null | tail -n +2 | grep -E "boltz2_|vina_|diffdock_|md_pbsa_" || echo "  None"

    echo ""
    echo "Pending jobs:"
    squeue -u $USER -t PD -o "  %.50j %.8T" 2>/dev/null | tail -n +2 | grep -E "boltz2_|vina_|diffdock_|md_pbsa_" || echo "  None"

    echo ""
}

show_help() {
    cat <<EOF
Master Automation Script for Screening Workflow
================================================

Usage: $0 [OPTIONS] [COMMAND]

Commands:
  all             Run all enabled workflows (default)
                  - MD+PBSA will wait for DiffDock to complete if WAIT_FOR_DIFFDOCK=1
  boltz2          Run only Boltz2 workflow
  vina            Run only Vina workflow
  diffdock        Run only DiffDock workflow
  md_pbsa         Run only MD+PBSA workflow
  status          Show current status of all workflows
  help            Show this help message

Options:
  --task-root DIR     Set task root directory (default: $TASK_ROOT)
  --proteins LIST     Set protein list (space-separated, default: $PROTEINS)
  --skip-boltz2       Skip Boltz2 workflow
  --skip-vina         Skip Vina workflow
  --skip-diffdock     Skip DiffDock workflow
  --skip-md-pbsa      Skip MD+PBSA workflow

Dependency Control:
  WAIT_FOR_DIFFDOCK=$WAIT_FOR_DIFFDOCK
    - When 1: MD+PBSA waits for DiffDock to complete before starting
    - When 0: MD+PBSA runs immediately without checking DiffDock

  WAIT_POLL_INTERVAL=$WAIT_POLL_INTERVAL seconds
    - How often to check DiffDock completion status

Examples:
  # Run all workflows with default settings (MD+PBSA waits for DiffDock)
  $0

  # Run only Boltz2 and Vina
  $0 --skip-diffdock --skip-md-pbsa

  # Run with custom task root and proteins
  $0 --task-root /path/to/data --proteins "PROTEIN1 PROTEIN2"

  # Run only DiffDock workflow
  $0 diffdock

  # Check status
  $0 status

Configuration:
  Task Root: $TASK_ROOT
  Proteins: $PROTEINS

  Workflows enabled:
    Boltz2: $([ $RUN_BOLTZ2 -eq 1 ] && echo "Yes" || echo "No")
    Vina: $([ $RUN_VINA -eq 1 ] && echo "Yes" || echo "No")
    DiffDock: $([ $RUN_DIFFDOCK -eq 1 ] && echo "Yes" || echo "No")
    MD+PBSA: $([ $RUN_MD_PBSA -eq 1 ] && echo "Yes" || echo "No")

  Dependency Control:
    Wait for DiffDock before MD+PBSA: $([ $WAIT_FOR_DIFFDOCK -eq 1 ] && echo "Yes" || echo "No")
    Check interval: ${WAIT_POLL_INTERVAL} seconds

Workflow Dependencies:
  - Boltz2: Independent (can run anytime)
  - Vina: Independent (can run anytime)
  - DiffDock: Requires Vina input preparation (split CSV files)
  - MD+PBSA: Requires DiffDock output (rank1.sdf files)

Notes:
  - When running 'all' command with WAIT_FOR_DIFFDOCK=1, the script will:
    1. Submit Boltz2, Vina, and DiffDock jobs
    2. Wait for all DiffDock jobs to complete
    3. Submit MD+PBSA jobs
  - Press Ctrl+C while waiting to skip MD+PBSA and exit
  - Individual workflow commands (boltz2, vina, diffdock, md_pbsa) ignore WAIT_FOR_DIFFDOCK

EOF
}

# ============================================================================
# Main Script
# ============================================================================

# Parse command line arguments
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
        --skip-boltz2)
            RUN_BOLTZ2=0
            shift
            ;;
        --skip-vina)
            RUN_VINA=0
            shift
            ;;
        --skip-diffdock)
            RUN_DIFFDOCK=0
            shift
            ;;
        --skip-md-pbsa)
            RUN_MD_PBSA=0
            shift
            ;;
        all|boltz2|vina|diffdock|md_pbsa|status|help)
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

# Handle commands
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

# Show banner
echo "=========================================="
echo "  Screening Workflow Automation"
echo "=========================================="
echo ""
echo "Task Root: $TASK_ROOT"
echo "Proteins: $PROTEINS"
echo "Command: $COMMAND"
echo ""

# Export configuration
export_config

# Check prerequisites
if ! check_prerequisites; then
    log_error "Prerequisites check failed. Please fix the issues and try again."
    exit 1
fi

# Execute workflows based on command
case $COMMAND in
    all)
        log_info "Running all enabled workflows"

        if [ $RUN_BOLTZ2 -eq 1 ]; then
            run_boltz2
        else
            echo "Skipping Boltz2 (disabled)"
            echo ""
        fi

        if [ $RUN_VINA -eq 1 ]; then
            run_vina
        else
            echo "Skipping Vina (disabled)"
            echo ""
        fi

        if [ $RUN_DIFFDOCK -eq 1 ]; then
            run_diffdock
        else
            echo "Skipping DiffDock (disabled)"
            echo ""
        fi

        if [ $RUN_MD_PBSA -eq 1 ]; then
            # Check if we should wait for DiffDock
            if [ $WAIT_FOR_DIFFDOCK -eq 1 ] && [ $RUN_DIFFDOCK -eq 1 ]; then
                echo ""
                log_info "MD+PBSA has dependency on DiffDock"

                # Wait for DiffDock to complete
                if wait_for_diffdock; then
                    run_md_pbsa
                else
                    log_error "DiffDock wait was interrupted. MD+PBSA not started."
                fi
            elif [ $WAIT_FOR_DIFFDOCK -eq 1 ] && [ $RUN_DIFFDOCK -eq 0 ]; then
                echo ""
                log_info "Checking DiffDock completion before MD+PBSA"

                # DiffDock was not run in this session, check if it's already complete
                if check_diffdock_complete; then
                    run_md_pbsa
                else
                    log_error "DiffDock is not complete. Please complete DiffDock first or disable WAIT_FOR_DIFFDOCK."
                    echo "Skipping MD+PBSA"
                    echo ""
                fi
            else
                # No dependency check, run immediately
                run_md_pbsa
            fi
        else
            echo "Skipping MD+PBSA (disabled)"
            echo ""
        fi
        ;;
    boltz2)
        run_boltz2
        ;;
    vina)
        run_vina
        ;;
    diffdock)
        run_diffdock
        ;;
    md_pbsa)
        run_md_pbsa
        ;;
esac

# Show final status
log_info "Workflow submission completed!"
echo "Use the following commands to monitor progress:"
echo ""
echo "  Overall status:        $0 status"
echo "  Boltz2 monitoring:     ${SCRIPT_DIR}/monitor_boltz2_jobs.sh"
echo "  Vina monitoring:       ${SCRIPT_DIR}/monitor_vina_jobs.sh"
echo "  DiffDock monitoring:   ${SCRIPT_DIR}/monitor_diffdock_jobs.sh"
echo "  MD+PBSA monitoring:    ${SCRIPT_DIR}/monitor_md_pbsa_jobs.sh"
echo ""
echo "  SLURM queue:           squeue -u \$USER"
echo ""
