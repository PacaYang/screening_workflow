#!/bin/bash
#
# Pipeline Utilities Library
# Common functions used across all pipeline stages
#

# Source logger (only if not already sourced)
if [ -z "$(type -t log_info)" ]; then
    LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${LIB_DIR}/logger.sh"
fi

# Check if a script exists
check_script_exists() {
    local script_path="$1"
    local script_name="$2"

    if [ ! -f "$script_path" ]; then
        log_error "Script not found: ${script_name} at ${script_path}"
        return 1
    fi

    if [ ! -x "$script_path" ]; then
        log_warn "Script not executable: ${script_name}, making it executable"
        chmod +x "$script_path"
    fi

    log_debug "Script found: ${script_name}"
    return 0
}

# Run command or dry-run
run_or_dry() {
    local dry_run="$1"
    shift
    local cmd="$@"

    if [ "$dry_run" = "1" ]; then
        log_info "[DRY RUN] Would execute: ${cmd}"
        return 0
    else
        log_debug "Executing: ${cmd}"
        eval "$cmd"
        return $?
    fi
}

# Collect all job files from multiple proteins
collect_all_job_files() {
    local task_root="$1"
    local proteins="$2"
    local job_file_pattern="$3"  # e.g., "initial_screening/GraphDTA/job_ids.txt"

    local all_job_files=()
    for protein in $proteins; do
        local job_file="${task_root}/${protein}/${job_file_pattern}"
        if [ -f "$job_file" ]; then
            all_job_files+=("$job_file")
        fi
    done

    echo "${all_job_files[@]}"
}

# Validate stage inputs
validate_stage_inputs() {
    local stage="$1"
    local task_root="$2"
    local proteins="$3"
    shift 3
    local required_files=("$@")

    log_debug "Validating inputs for stage ${stage}"

    for protein in $proteins; do
        for file_pattern in "${required_files[@]}"; do
            local file_path="${task_root}/${protein}/${file_pattern}"
            if [ ! -e "$file_path" ]; then
                log_error "Required input missing: ${file_path}"
                return 1
            fi
        done
    done

    log_debug "All inputs validated for stage ${stage}"
    return 0
}

# Validate stage outputs
validate_stage_outputs() {
    local stage="$1"
    local task_root="$2"
    local proteins="$3"
    shift 3
    local required_files=("$@")

    log_debug "Validating outputs for stage ${stage}"

    for protein in $proteins; do
        for file_pattern in "${required_files[@]}"; do
            local file_path="${task_root}/${protein}/${file_pattern}"
            if [ ! -e "$file_path" ]; then
                log_error "Required output missing: ${file_path}"
                return 1
            fi
        done
    done

    log_debug "All outputs validated for stage ${stage}"
    return 0
}

# Check if stage is complete
is_stage_complete() {
    local task_root="$1"
    local proteins="$2"
    shift 2
    local completion_markers=("$@")

    for protein in $proteins; do
        for marker in "${completion_markers[@]}"; do
            local marker_path="${task_root}/${protein}/${marker}"
            if [ ! -e "$marker_path" ]; then
                return 1
            fi
        done
    done

    return 0
}

# Export environment variables for pipeline
export_pipeline_env() {
    local task_root="$1"
    local proteins="$2"

    export MASTER_TASK_ROOT="$task_root"
    export MASTER_PROTEINS="$proteins"

    log_debug "Exported environment: MASTER_TASK_ROOT=${task_root}, MASTER_PROTEINS=${proteins}"
}

# Parse boolean flag
parse_bool_flag() {
    local value="$1"
    case "$value" in
        1|true|yes|on) echo "1" ;;
        0|false|no|off) echo "0" ;;
        *) echo "0" ;;
    esac
}
