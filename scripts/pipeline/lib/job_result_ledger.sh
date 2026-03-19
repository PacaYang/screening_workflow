#!/bin/bash
#
# Job Result Ledger Utilities
# Append per-job success/fail records in a concurrency-safe way.
#

# Append one job result line to a TSV ledger file.
# Columns:
#   timestamp_utc \t job_id \t protein \t job_key \t status \t completed_count \t expected_count
# Args:
#   ledger_file job_id protein job_key status completed_count expected_count
append_job_result() {
    local ledger_file="$1"
    local job_id="$2"
    local protein="$3"
    local job_key="$4"
    local status="$5"
    local completed_count="$6"
    local expected_count="$7"

    mkdir -p "$(dirname "$ledger_file")"
    touch "$ledger_file"

    local timestamp
    timestamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    local line
    line="${timestamp}\t${job_id}\t${protein}\t${job_key}\t${status}\t${completed_count}\t${expected_count}"

    if command -v flock >/dev/null 2>&1; then
        local lock_file="${ledger_file}.lock"
        : > "$lock_file"
        exec 9>>"$lock_file"
        flock 9
        printf '%b\n' "$line" >> "$ledger_file"
        flock -u 9
        exec 9>&-
    else
        printf '%b\n' "$line" >> "$ledger_file"
    fi
}
