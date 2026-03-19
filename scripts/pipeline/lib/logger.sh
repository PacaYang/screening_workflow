#!/bin/bash
#
# Logging Library for Pipeline
# Provides consistent logging across all pipeline components
#

# Global log file (set by caller)
LOG_FILE="${LOG_FILE:-}"
DEBUG="${DEBUG:-0}"

# Color codes for terminal output
COLOR_RESET="\033[0m"
COLOR_RED="\033[31m"
COLOR_YELLOW="\033[33m"
COLOR_GREEN="\033[32m"
COLOR_BLUE="\033[34m"
COLOR_GRAY="\033[90m"

# Log levels
LOG_LEVEL_DEBUG=0
LOG_LEVEL_INFO=1
LOG_LEVEL_WARN=2
LOG_LEVEL_ERROR=3

# Get timestamp
get_timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

# Write to log file and stdout
write_log() {
    local level="$1"
    local message="$2"
    local color="$3"

    local timestamp
    timestamp=$(get_timestamp)
    local log_line="[${timestamp}] [${level}] ${message}"

    # Write to log file if configured
    if [ -n "$LOG_FILE" ]; then
        echo "$log_line" >> "$LOG_FILE"
    fi

    # Write to stdout with color
    if [ -t 1 ]; then
        echo -e "${color}[${timestamp}] [${level}]${COLOR_RESET} ${message}"
    else
        echo "$log_line"
    fi
}

# Log debug message (only if DEBUG=1)
log_debug() {
    if [ "$DEBUG" = "1" ]; then
        write_log "DEBUG" "$1" "$COLOR_GRAY"
    fi
}

# Log info message
log_info() {
    write_log "INFO" "$1" "$COLOR_BLUE"
}

# Log warning message
log_warn() {
    write_log "WARN" "$1" "$COLOR_YELLOW"
}

# Log error message
log_error() {
    write_log "ERROR" "$1" "$COLOR_RED"
}

# Log stage start
log_stage_start() {
    local stage_name="$1"
    write_log "STAGE" "Starting: ${stage_name}" "$COLOR_GREEN"
}

# Log stage complete
log_stage_complete() {
    local stage_name="$1"
    write_log "STAGE" "Completed: ${stage_name}" "$COLOR_GREEN"
}

# Log job submission
log_job_submit() {
    local job_id="$1"
    local job_name="$2"
    write_log "JOB" "Submitted: ${job_name} (ID: ${job_id})" "$COLOR_BLUE"
}

# Log job completion
log_job_complete() {
    local job_id="$1"
    local job_name="$2"
    local status="$3"
    if [ "$status" = "0" ] || [ "$status" = "COMPLETED" ]; then
        write_log "JOB" "Completed: ${job_name} (ID: ${job_id})" "$COLOR_GREEN"
    else
        write_log "JOB" "Failed: ${job_name} (ID: ${job_id}, Status: ${status})" "$COLOR_RED"
    fi
}
