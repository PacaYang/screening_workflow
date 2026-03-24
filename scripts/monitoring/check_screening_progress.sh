#!/bin/bash
#
# Comprehensive Screening Progress Checker
# Checks completion and failure status for AF3, Vina, DiffDock, and MD+PBSA
#

set -e

# ============================================================================
# Parse command line arguments
# ============================================================================

show_usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Options:
  --task-root DIR       Set task root directory
  --proteins LIST       Set protein list (space-separated, quoted)
  -h, --help           Show this help message

Examples:
  # Use environment variables (MASTER_TASK_ROOT, MASTER_PROTEINS)
  $0

  # Specify task root and one protein
  $0 --task-root /home/ubuntu/snake_test --proteins "JAK1JH1"

  # Specify multiple proteins
  $0 --task-root /path/to/data --proteins "JAK1 JAK2 JAK3"

Environment Variables (used as defaults):
  MASTER_TASK_ROOT     Task root directory (default: /home/ubuntu/snake_test)
  MASTER_PROTEINS      Space-separated protein list (default: JAK1JH1)

EOF
}

# Default values from environment or hardcoded defaults
TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"
PROTEINS="${MASTER_PROTEINS:-JAK1JH1}"

# Parse command line arguments
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
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run '$0 --help' for usage information"
            exit 1
            ;;
    esac
done

# ============================================================================
# Functions
# ============================================================================

# Color codes for better readability
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_header() {
    echo ""
    echo "=========================================="
    echo "$1"
    echo "=========================================="
}

print_section() {
    echo ""
    echo "------------------------------------------"
    echo "$1"
    echo "------------------------------------------"
}

# Check AF3 progress
check_af3() {
    local protein=$1
    local n_compounds=$2

    print_section "AlphaFold 3"

    local af3_dir="${TASK_ROOT}/${protein}/fine_screening/AF3/output"

    if [ ! -d "$af3_dir" ]; then
        echo "  Status: Not started (output directory not found)"
        return
    fi

    local completed=0
    local missing=0

    # AF3 uses lowercase protein names in output
    local protein_lower=$(echo "$protein" | tr '[:upper:]' '[:lower:]')

    for i in $(seq 0 $((n_compounds - 1))); do
        local compound_name="${protein_lower}${i}"
        local json_file="${af3_dir}/${compound_name}/${compound_name}_summary_confidences.json"

        if [ -f "$json_file" ]; then
            completed=$((completed + 1))
        else
            missing=$((missing + 1))
        fi
    done

    local success_rate=0
    if [ $n_compounds -gt 0 ]; then
        success_rate=$((completed * 100 / n_compounds))
    fi

    echo "  Total compounds: $n_compounds"
    echo -e "  ${GREEN}Completed: $completed${NC}"
    echo -e "  ${YELLOW}Pending: $missing${NC}"
    echo "  Success rate: ${success_rate}%"

    # Show sample successful output
    if [ $completed -gt 0 ]; then
        local sample_json=$(find "$af3_dir" -name "*_summary_confidences.json" -type f | head -1)
        if [ -f "$sample_json" ]; then
            echo ""
            echo "  Sample output: $(basename $(dirname "$sample_json"))"
        fi
    fi
}

# Check Boltz2 progress
check_boltz2() {
    local protein=$1
    local n_compounds=$2

    print_section "Boltz2"

    local boltz2_dir="${TASK_ROOT}/${protein}/fine_screening/Boltz2/output"
    local token_dir="${boltz2_dir}/token"

    if [ ! -d "$boltz2_dir" ]; then
        echo "  Status: Not started (output directory not found)"
        return
    fi

    # Boltz2 uses 40 batches by default
    local n_batches=40
    local completed_batches=0
    local missing_batches=0

    # Check batch completion tokens
    if [ -d "$token_dir" ]; then
        completed_batches=$(find "$token_dir" -name "batch_*.done" -type f 2>/dev/null | wc -l)
    fi
    missing_batches=$((n_batches - completed_batches))

    # Check for compressed tar files
    local tar_files=0
    if [ -d "$boltz2_dir" ]; then
        tar_files=$(find "$boltz2_dir" -maxdepth 1 -name "batch_*.tar.gz" -type f 2>/dev/null | wc -l)
    fi

    local batch_success_rate=0
    if [ $n_batches -gt 0 ]; then
        batch_success_rate=$((completed_batches * 100 / n_batches))
    fi

    # Estimate compounds per batch
    local compounds_per_batch=$((n_compounds / n_batches))
    if [ $compounds_per_batch -eq 0 ]; then
        compounds_per_batch=1
    fi
    local estimated_completed=$((completed_batches * compounds_per_batch))

    echo "  Total compounds: $n_compounds"
    echo "  Batch processing: $n_batches batches"
    echo ""
    echo "  Batch Status:"
    echo -e "    ${GREEN}Completed batches: $completed_batches${NC}"
    echo -e "    ${YELLOW}Pending batches: $missing_batches${NC}"
    echo "    Batch completion rate: ${batch_success_rate}%"
    echo ""
    echo "  Output Files:"
    echo "    Compressed archives (tar.gz): $tar_files"
    echo ""
    echo "  Estimated Progress:"
    echo "    ~$estimated_completed / $n_compounds compounds completed"
    echo "    (~$((completed_batches * 100 / n_batches))% complete)"

    # Show sample output
    if [ $completed_batches -gt 0 ]; then
        local sample_token=$(find "$token_dir" -name "batch_*.done" -type f 2>/dev/null | head -1)
        if [ -f "$sample_token" ]; then
            echo ""
            echo "  Sample completed batch: $(basename "$sample_token" .done)"
        fi

        local sample_tar=$(find "$boltz2_dir" -maxdepth 1 -name "batch_*.tar.gz" -type f 2>/dev/null | head -1)
        if [ -f "$sample_tar" ]; then
            echo "  Sample archive: $(basename "$sample_tar") ($(du -h "$sample_tar" | cut -f1))"
        fi
    fi

    # Check SLURM jobs for this protein's Boltz2 batches
    local running_jobs=$(squeue -u $USER --states=R 2>/dev/null | grep "boltz2_${protein}_" | wc -l)
    local pending_jobs=$(squeue -u $USER --states=PD 2>/dev/null | grep "boltz2_${protein}_" | wc -l)

    if [ $running_jobs -gt 0 ] || [ $pending_jobs -gt 0 ]; then
        echo ""
        echo "  Active SLURM Jobs:"
        echo "    Running: $running_jobs"
        echo "    Pending: $pending_jobs"
    fi
}

# Check Vina progress
check_vina() {
    local protein=$1
    local n_compounds=$2

    print_section "AutoDock Vina"

    local vina_output="${TASK_ROOT}/${protein}/fine_screening/Vina/output"
    local vina_input="${TASK_ROOT}/${protein}/fine_screening/Vina/input"

    if [ ! -d "$vina_output" ]; then
        echo "  Status: Not started (output directory not found)"
        return
    fi

    if [ ! -d "$vina_input" ]; then
        echo "  Status: Cannot determine progress (input directory not found)"
        return
    fi

    # Count parts
    local n_parts=$(ls "$vina_input"/*.csv 2>/dev/null | wc -l)

    if [ $n_parts -eq 0 ]; then
        echo "  Status: No input parts found"
        return
    fi

    local completed=0
    local missing=0
    local total_ligands=0

    # Check each part
    for part_csv in "$vina_input"/*.csv; do
        local part_name=$(basename "$part_csv" .csv)
        local part_dir="${vina_output}/${part_name}"

        # Count expected ligands in this part
        local n_ligands=$(($(wc -l < "$part_csv") - 1))  # Subtract header
        total_ligands=$((total_ligands + n_ligands))

        if [ ! -d "$part_dir" ]; then
            missing=$((missing + n_ligands))
            continue
        fi

        # Check each ligand in the part
        for j in $(seq 0 $((n_ligands - 1))); do
            local affinity_file="${part_dir}/lig${j}/docking_affinities.txt"

            if [ -f "$affinity_file" ]; then
                completed=$((completed + 1))
            else
                missing=$((missing + 1))
            fi
        done
    done

    local success_rate=0
    if [ $total_ligands -gt 0 ]; then
        success_rate=$((completed * 100 / total_ligands))
    fi

    echo "  Total compounds: $total_ligands (across $n_parts parts)"
    echo -e "  ${GREEN}Completed: $completed${NC}"
    echo -e "  ${YELLOW}Pending: $missing${NC}"
    echo "  Success rate: ${success_rate}%"

    # Show sample output
    if [ $completed -gt 0 ]; then
        local sample_file=$(find "$vina_output" -name "docking_affinities.txt" -type f | head -1)
        if [ -f "$sample_file" ]; then
            echo ""
            echo "  Sample output: $(dirname "$sample_file" | sed "s|$vina_output/||")"
        fi
    fi
}

# Check DiffDock progress
check_diffdock() {
    local protein=$1
    local n_compounds=$2

    print_section "DiffDock"

    local diffdock_dir="${TASK_ROOT}/${protein}/fine_screening/PBSA/DiffDock/output"

    if [ ! -d "$diffdock_dir" ]; then
        echo "  Status: Not started (output directory not found)"
        return
    fi

    local completed=0
    local missing=0

    for i in $(seq 0 $((n_compounds - 1))); do
        local compound_dir="${diffdock_dir}/${protein}_${i}"
        local sdf_file="${compound_dir}/rank1.sdf"

        if [ -f "$sdf_file" ]; then
            completed=$((completed + 1))
        else
            missing=$((missing + 1))
        fi
    done

    local success_rate=0
    if [ $n_compounds -gt 0 ]; then
        success_rate=$((completed * 100 / n_compounds))
    fi

    echo "  Total compounds: $n_compounds"
    echo -e "  ${GREEN}Completed: $completed${NC}"
    echo -e "  ${YELLOW}Pending: $missing${NC}"
    echo "  Success rate: ${success_rate}%"

    # Show sample output
    if [ $completed -gt 0 ]; then
        local sample_sdf=$(find "$diffdock_dir" -name "rank1.sdf" -type f | head -1)
        if [ -f "$sample_sdf" ]; then
            local sample_dir=$(dirname "$sample_sdf")
            echo ""
            echo "  Sample output: $(basename "$sample_dir")"
            echo "  File size: $(du -h "$sample_sdf" | cut -f1)"
        fi
    fi
}

# Check MD+PBSA progress
check_md_pbsa() {
    local protein=$1
    local n_compounds=$2

    print_section "MD + PBSA"

    local md_dir="${TASK_ROOT}/${protein}/fine_screening/PBSA/PBSA/MD"
    local pbsa_dir="${TASK_ROOT}/${protein}/fine_screening/PBSA/PBSA/PBSA"

    # Check MD phase
    echo "  MD Phase:"
    if [ ! -d "$md_dir" ]; then
        echo "    Status: Not started (MD directory not found)"
        local md_completed=0
        local md_missing=$n_compounds
    else
        local md_completed=0
        local md_missing=0

        for i in $(seq 0 $((n_compounds - 1))); do
            local md_compound_dir="${md_dir}/${protein}_${i}"
            local trajectory="${md_compound_dir}/T298.xtc"

            if [ -f "$trajectory" ]; then
                md_completed=$((md_completed + 1))
            else
                md_missing=$((md_missing + 1))
            fi
        done
    fi

    local md_success_rate=0
    if [ $n_compounds -gt 0 ]; then
        md_success_rate=$((md_completed * 100 / n_compounds))
    fi

    echo -e "    ${GREEN}Completed: $md_completed${NC}"
    echo -e "    ${YELLOW}Pending: $md_missing${NC}"
    echo "    Success rate: ${md_success_rate}%"

    # Check PBSA phase
    echo ""
    echo "  PBSA Phase:"
    if [ ! -d "$pbsa_dir" ]; then
        echo "    Status: Not started (PBSA directory not found)"
        local pbsa_completed=0
        local pbsa_missing=$n_compounds
    else
        local pbsa_completed=0
        local pbsa_missing=0

        for i in $(seq 0 $((n_compounds - 1))); do
            local pbsa_compound_dir="${pbsa_dir}/${protein}_${i}"
            local pbsa_dat="${pbsa_compound_dir}/FINAL_RESULTS_MMPBSA.dat"

            if [ -f "$pbsa_dat" ]; then
                pbsa_completed=$((pbsa_completed + 1))
            else
                pbsa_missing=$((pbsa_missing + 1))
            fi
        done
    fi

    local pbsa_success_rate=0
    if [ $n_compounds -gt 0 ]; then
        pbsa_success_rate=$((pbsa_completed * 100 / n_compounds))
    fi

    echo -e "    ${GREEN}Completed: $pbsa_completed${NC}"
    echo -e "    ${YELLOW}Pending: $pbsa_missing${NC}"
    echo "    Success rate: ${pbsa_success_rate}%"

    # Overall summary
    echo ""
    echo "  Overall (MD → PBSA pipeline):"
    echo "    Total compounds: $n_compounds"
    echo "    Fully completed (MD+PBSA): $pbsa_completed"
    echo "    In progress (MD done, PBSA pending): $((md_completed - pbsa_completed))"

    # Show sample output
    if [ $pbsa_completed -gt 0 ]; then
        local sample_dat=$(find "$pbsa_dir" -name "FINAL_RESULTS_MMPBSA.dat" -type f | head -1)
        if [ -f "$sample_dat" ]; then
            echo ""
            echo "  Sample PBSA output: $(basename $(dirname "$sample_dat"))"
        fi
    fi
}

# Check SLURM jobs
check_slurm_jobs() {
    local protein=$1

    print_section "SLURM Job Status"

    local af3_running=$(squeue -u $USER --states=R -o "%.100j" 2>/dev/null | grep "af3_${protein}_" | wc -l)
    local af3_pending=$(squeue -u $USER --states=PD -o "%.100j" 2>/dev/null | grep "af3_${protein}_" | wc -l)

    local vina_running=$(squeue -u $USER --states=R -o "%.100j" 2>/dev/null | grep "vina_${protein}_" | wc -l)
    local vina_pending=$(squeue -u $USER --states=PD -o "%.100j" 2>/dev/null | grep "vina_${protein}_" | wc -l)

    local diffdock_running=$(squeue -u $USER --states=R -o "%.100j" 2>/dev/null | grep "diffdock_${protein}_" | wc -l)
    local diffdock_pending=$(squeue -u $USER --states=PD -o "%.100j" 2>/dev/null | grep "diffdock_${protein}_" | wc -l)

    local md_pbsa_running=$(squeue -u $USER --states=R -o "%.100j" 2>/dev/null | grep "md_pbsa_${protein}_" | wc -l)
    local md_pbsa_pending=$(squeue -u $USER --states=PD -o "%.100j" 2>/dev/null | grep "md_pbsa_${protein}_" | wc -l)

    echo "  AF3:       Running: $af3_running, Pending: $af3_pending"
    echo "  Vina:      Running: $vina_running, Pending: $vina_pending"
    echo "  DiffDock:  Running: $diffdock_running, Pending: $diffdock_pending"
    echo "  MD+PBSA:   Running: $md_pbsa_running, Pending: $md_pbsa_pending"

    local total_running=$((af3_running + vina_running + diffdock_running + md_pbsa_running))
    local total_pending=$((af3_pending + vina_pending + diffdock_pending + md_pbsa_pending))

    echo ""
    echo "  Total:     Running: $total_running, Pending: $total_pending"
}

# Generate summary report
generate_summary() {
    print_header "SUMMARY REPORT"

    local total_proteins=0
    local total_compounds=0

    for protein in $PROTEINS; do
        local selected_csv="${TASK_ROOT}/${protein}/initial_screening/selected.csv"
        if [ -f "$selected_csv" ]; then
            local n=$(tail -n +2 "$selected_csv" | wc -l)
            total_compounds=$((total_compounds + n))
            total_proteins=$((total_proteins + 1))
        fi
    done

    echo "  Total proteins: $total_proteins"
    echo "  Total compounds: $total_compounds"
    echo ""
    echo "  Task root: $TASK_ROOT"
    echo "  Proteins: $PROTEINS"
}

# ============================================================================
# Main Script
# ============================================================================

echo "=========================================="
echo "  Screening Progress Report"
echo "=========================================="
echo "  Generated at: $(date)"
echo "=========================================="

# Check if task root exists
if [ ! -d "$TASK_ROOT" ]; then
    echo ""
    echo "ERROR: Task root directory not found: $TASK_ROOT"
    echo "Please check the path and try again."
    exit 1
fi

# Generate summary
generate_summary

# Check each protein
for protein in $PROTEINS; do
    print_header "Protein: $protein"

    # Get number of compounds
    selected_csv="${TASK_ROOT}/${protein}/initial_screening/selected.csv"

    if [ ! -f "$selected_csv" ]; then
        echo "  ERROR: Selected compounds file not found: $selected_csv"
        echo "  Skipping this protein."
        continue
    fi

    n_compounds=$(tail -n +2 "$selected_csv" | wc -l)
    echo "  Total compounds in screening: $n_compounds"

    # Check each algorithm
    check_boltz2 "$protein" "$n_compounds"
    check_af3 "$protein" "$n_compounds"
    check_vina "$protein" "$n_compounds"
    check_diffdock "$protein" "$n_compounds"
    check_md_pbsa "$protein" "$n_compounds"
    check_slurm_jobs "$protein"
done

# Final summary
print_header "Screening Check Complete"
echo "  Timestamp: $(date)"
echo "=========================================="
echo ""
echo "For detailed monitoring of individual workflows, use:"
echo "  - Boltz2: ./monitor_boltz2_jobs.sh"
echo "  - AF3: ./monitor_af3_jobs.sh"
echo "  - Vina: ./monitor_vina_jobs.sh"
echo "  - DiffDock: ./monitor_diffdock_jobs.sh"
echo "  - MD+PBSA: ./monitor_md_pbsa_jobs.sh"
echo ""
