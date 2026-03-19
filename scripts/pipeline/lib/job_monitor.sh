#!/bin/bash
#
# Job Monitoring Library for Pipeline
# SLURM job submission and monitoring functions
#

# Source dependencies (only if not already sourced)
if [ -z "$(type -t log_info)" ]; then
    LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${LIB_DIR}/logger.sh"
    source "${LIB_DIR}/pipeline_utils.sh"
fi

# Wait for SLURM jobs to complete
# Args: poll_interval job_file1 [job_file2 ...]
wait_for_slurm_jobs() {
    local poll_interval="$1"
    shift
    local job_files=("$@")

    if [ ${#job_files[@]} -eq 0 ]; then
        log_warn "No job files provided to wait_for_slurm_jobs"
        return 0
    fi

    log_info "Waiting for SLURM jobs to complete (polling every ${poll_interval}s)"
    log_debug "Monitoring ${#job_files[@]} job files"

    local all_done=0
    while [ $all_done -eq 0 ]; do
        all_done=1

        for job_file in "${job_files[@]}"; do
            if [ ! -f "$job_file" ]; then
                log_warn "Job file not found: ${job_file}"
                continue
            fi

            while IFS= read -r job_id; do
                [ -z "$job_id" ] && continue

                if squeue -j "$job_id" &>/dev/null; then
                    all_done=0
                    log_debug "Job ${job_id} still running"
                    break
                fi
            done < "$job_file"

            if [ $all_done -eq 0 ]; then
                break
            fi
        done

        if [ $all_done -eq 0 ]; then
            sleep "$poll_interval"
        fi
    done

    log_info "All SLURM jobs completed"
}

# Submit job and track ID
submit_job_with_tracking() {
    local job_script="$1"
    local job_name="$2"
    local output_file="$3"

    if [ ! -f "$job_script" ]; then
        log_error "Job script not found: ${job_script}"
        return 1
    fi

    local job_id
    job_id=$(sbatch "$job_script" | awk '{print $NF}')

    if [ -z "$job_id" ]; then
        log_error "Failed to submit job: ${job_name}"
        return 1
    fi

    log_job_submit "$job_id" "$job_name"

    # Save job ID to file
    if [ -n "$output_file" ]; then
        echo "$job_id" >> "$output_file"
    fi

    echo "$job_id"
}

# Get job status
get_job_status() {
    local job_id="$1"

    if ! squeue -j "$job_id" &>/dev/null; then
        # Job not in queue, check sacct for completion status
        local status
        status=$(sacct -j "$job_id" --format=State --noheader | head -1 | tr -d ' ')
        echo "${status:-COMPLETED}"
    else
        echo "RUNNING"
    fi
}

# Get failed jobs from a list
get_failed_jobs() {
    local job_file="$1"

    if [ ! -f "$job_file" ]; then
        return 0
    fi

    local failed_jobs=()
    while IFS= read -r job_id; do
        [ -z "$job_id" ] && continue

        local status
        status=$(get_job_status "$job_id")

        if [[ "$status" != "COMPLETED" && "$status" != "RUNNING" ]]; then
            failed_jobs+=("$job_id")
        fi
    done < "$job_file"

    echo "${failed_jobs[@]}"
}

# Cancel jobs
cancel_jobs() {
    local job_file="$1"

    if [ ! -f "$job_file" ]; then
        log_warn "Job file not found: ${job_file}"
        return 0
    fi

    log_info "Cancelling jobs from ${job_file}"

    while IFS= read -r job_id; do
        [ -z "$job_id" ] && continue

        if squeue -j "$job_id" &>/dev/null; then
            scancel "$job_id"
            log_info "Cancelled job: ${job_id}"
        fi
    done < "$job_file"
}

# Wait for jobs with progress reporting
wait_for_jobs_with_progress() {
    local poll_interval="$1"
    shift
    local job_files=("$@")

    if [ ${#job_files[@]} -eq 0 ]; then
        return 0
    fi

    # Count total jobs
    local total_jobs=0
    for job_file in "${job_files[@]}"; do
        if [ -f "$job_file" ]; then
            total_jobs=$((total_jobs + $(wc -l < "$job_file")))
        fi
    done

    log_info "Monitoring ${total_jobs} jobs (polling every ${poll_interval}s)"

    local all_done=0
    while [ $all_done -eq 0 ]; do
        all_done=1
        local running_jobs=0

        for job_file in "${job_files[@]}"; do
            if [ ! -f "$job_file" ]; then
                continue
            fi

            while IFS= read -r job_id; do
                [ -z "$job_id" ] && continue

                if squeue -j "$job_id" &>/dev/null; then
                    all_done=0
                    running_jobs=$((running_jobs + 1))
                fi
            done < "$job_file"
        done

        if [ $all_done -eq 0 ]; then
            local completed=$((total_jobs - running_jobs))
            local progress=$((completed * 100 / total_jobs))
            log_info "Progress: ${completed}/${total_jobs} jobs completed (${progress}%)"
            sleep "$poll_interval"
        fi
    done

    log_info "All jobs completed"
}
