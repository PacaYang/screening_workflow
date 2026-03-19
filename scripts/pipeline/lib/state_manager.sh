#!/bin/bash
#
# State Management Library for Pipeline
# Tracks pipeline progress using JSON state files
#

# Source logger (only if not already sourced)
if [ -z "$(type -t log_info)" ]; then
    LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${LIB_DIR}/logger.sh"
fi

# Global state variables
PIPELINE_ID=""
STATE_FILE=""
STATE_DIR=""

# Initialize pipeline state
# Creates a new state file with unique pipeline_id
init_pipeline_state() {
    local task_root="$1"
    local proteins="$2"

    # Create state directory
    STATE_DIR="${task_root}/.pipeline_state"
    mkdir -p "$STATE_DIR"

    # Generate unique pipeline ID
    PIPELINE_ID="$(date +%Y%m%d_%H%M%S)_$(echo "$proteins" | tr ' ' '_' | head -c 20)"
    STATE_FILE="${STATE_DIR}/${PIPELINE_ID}.json"

    # Create initial state
    cat > "$STATE_FILE" <<EOF
{
  "pipeline_id": "${PIPELINE_ID}",
  "task_root": "${task_root}",
  "proteins": $(echo "$proteins" | jq -R 'split(" ")'),
  "start_time": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "end_time": null,
  "status": "running",
  "stages": {
    "1": {"status": "pending", "start": null, "end": null, "exit_code": null, "job_ids": []},
    "2": {"status": "pending", "start": null, "end": null, "exit_code": null, "job_ids": []},
    "3": {"status": "pending", "start": null, "end": null, "exit_code": null, "job_ids": []},
    "4": {"status": "pending", "start": null, "end": null, "exit_code": null, "job_ids": []},
    "5": {"status": "pending", "start": null, "end": null, "exit_code": null, "job_ids": []},
    "6": {"status": "pending", "start": null, "end": null, "exit_code": null, "job_ids": []}
  },
  "config": {}
}
EOF

    log_info "Initialized pipeline state: ${PIPELINE_ID}"
    log_debug "State file: ${STATE_FILE}"
}

# Load existing pipeline state
load_pipeline_state() {
    local pipeline_id="$1"
    local task_root="$2"

    STATE_DIR="${task_root}/.pipeline_state"
    STATE_FILE="${STATE_DIR}/${pipeline_id}.json"

    if [ ! -f "$STATE_FILE" ]; then
        log_error "State file not found: ${STATE_FILE}"
        return 1
    fi

    PIPELINE_ID="$pipeline_id"
    log_info "Loaded pipeline state: ${PIPELINE_ID}"
}

# Update stage status
update_stage_status() {
    local stage="$1"
    local status="$2"  # pending, running, completed, failed
    local exit_code="${3:-}"

    if [ ! -f "$STATE_FILE" ]; then
        log_error "State file not found: ${STATE_FILE}"
        return 1
    fi

    local timestamp
    timestamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)

    # Update status and timestamps
    local tmp_file="${STATE_FILE}.tmp"
    if [ "$status" = "running" ]; then
        jq ".stages.\"${stage}\".status = \"running\" | .stages.\"${stage}\".start = \"${timestamp}\"" "$STATE_FILE" > "$tmp_file"
    elif [ "$status" = "completed" ] || [ "$status" = "failed" ]; then
        jq ".stages.\"${stage}\".status = \"${status}\" | .stages.\"${stage}\".end = \"${timestamp}\" | .stages.\"${stage}\".exit_code = ${exit_code:-0}" "$STATE_FILE" > "$tmp_file"
    else
        jq ".stages.\"${stage}\".status = \"${status}\"" "$STATE_FILE" > "$tmp_file"
    fi

    mv "$tmp_file" "$STATE_FILE"
    log_debug "Updated stage ${stage} status: ${status}"
}

# Get stage status
get_stage_status() {
    local stage="$1"

    if [ ! -f "$STATE_FILE" ]; then
        echo "unknown"
        return 1
    fi

    jq -r ".stages.\"${stage}\".status" "$STATE_FILE"
}

# Record job IDs for a stage
record_job_ids() {
    local stage="$1"
    shift
    local job_ids=("$@")

    if [ ! -f "$STATE_FILE" ]; then
        log_error "State file not found: ${STATE_FILE}"
        return 1
    fi

    # Convert array to JSON array
    local json_array
    json_array=$(printf '%s\n' "${job_ids[@]}" | jq -R . | jq -s .)

    local tmp_file="${STATE_FILE}.tmp"
    jq ".stages.\"${stage}\".job_ids = ${json_array}" "$STATE_FILE" > "$tmp_file"
    mv "$tmp_file" "$STATE_FILE"

    log_debug "Recorded ${#job_ids[@]} job IDs for stage ${stage}"
}

# Get job IDs for a stage
get_job_ids() {
    local stage="$1"

    if [ ! -f "$STATE_FILE" ]; then
        return 1
    fi

    jq -r ".stages.\"${stage}\".job_ids[]" "$STATE_FILE" 2>/dev/null
}

# Get pipeline progress (percentage)
get_pipeline_progress() {
    if [ ! -f "$STATE_FILE" ]; then
        echo "0"
        return 1
    fi

    local completed
    completed=$(jq '[.stages[] | select(.status == "completed")] | length' "$STATE_FILE")
    local total=6
    echo $((completed * 100 / total))
}

# Mark pipeline as complete
mark_pipeline_complete() {
    local status="$1"  # completed or failed

    if [ ! -f "$STATE_FILE" ]; then
        return 1
    fi

    local timestamp
    timestamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)

    local tmp_file="${STATE_FILE}.tmp"
    jq ".status = \"${status}\" | .end_time = \"${timestamp}\"" "$STATE_FILE" > "$tmp_file"
    mv "$tmp_file" "$STATE_FILE"

    log_info "Pipeline ${status}: ${PIPELINE_ID}"
}

# Save configuration to state
save_config_to_state() {
    local config_json="$1"

    if [ ! -f "$STATE_FILE" ]; then
        return 1
    fi

    local tmp_file="${STATE_FILE}.tmp"
    jq ".config = ${config_json}" "$STATE_FILE" > "$tmp_file"
    mv "$tmp_file" "$STATE_FILE"
}

# Get latest pipeline ID for a task root
get_latest_pipeline_id() {
    local task_root="$1"
    local state_dir="${task_root}/.pipeline_state"

    if [ ! -d "$state_dir" ]; then
        return 1
    fi

    ls -t "${state_dir}"/*.json 2>/dev/null | head -1 | xargs basename | sed 's/\.json$//'
}
