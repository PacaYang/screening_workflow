#!/bin/bash
#
# Pipeline Progress Monitor
# Writes a live text status file from method ledgers and job lists.
#

# Source logger (only if not already sourced)
if [ -z "$(type -t log_info)" ]; then
    LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${LIB_DIR}/logger.sh"
fi

PROGRESS_MONITOR_PID=""
PROGRESS_MONITOR_TASK_ROOT=""
PROGRESS_MONITOR_PROTEINS=""
PROGRESS_MONITOR_PIPELINE_TYPE=""
PROGRESS_MONITOR_INTERVAL=300
PROGRESS_STATUS_FILE=""

count_nonempty_lines() {
    local file="$1"
    if [ ! -f "$file" ]; then
        echo 0
        return 0
    fi
    local n
    n=$(grep -cve '^[[:space:]]*$' "$file" 2>/dev/null || true)
    echo "${n:-0}"
}

count_input_chunks() {
    local task_root="$1"
    local proteins="$2"
    local total=0

    for protein in $proteins; do
        local input_dir="${task_root}/${protein}/initial_screening/inputs"
        if [ -d "$input_dir" ]; then
            local n
            n=$(find "$input_dir" -maxdepth 1 -type f -name "input_*.csv" 2>/dev/null | wc -l)
            total=$((total + n))
        fi
    done

    echo "$total"
}

count_initial_predictions() {
    local task_root="$1"
    local proteins="$2"
    local method="$3"
    local total=0

    for protein in $proteins; do
        local method_dir="${task_root}/${protein}/initial_screening/${method}"
        if [ -d "$method_dir" ]; then
            local n
            n=$(find "$method_dir" -maxdepth 1 -type f -name "prediction_*.csv" 2>/dev/null | wc -l)
            total=$((total + n))
        fi
    done

    echo "$total"
}

format_simple_status() {
    local completed="$1"
    local expected="$2"

    if [ "$expected" -eq 0 ]; then
        echo "Pending"
    elif [ "$completed" -ge "$expected" ]; then
        echo "Done"
    elif [ "$completed" -eq 0 ]; then
        echo "Pending"
    else
        echo "Processing ${completed}/${expected}"
    fi
}

collect_method_expected_jobs() {
    local task_root="$1"
    local proteins="$2"
    local job_file_rel="$3"
    local expected=0

    for protein in $proteins; do
        local job_file="${task_root}/${protein}/${job_file_rel}"
        expected=$((expected + $(count_nonempty_lines "$job_file")))
    done

    echo "$expected"
}

collect_method_ledger_counts() {
    local task_root="$1"
    local proteins="$2"
    local ledger_rel="$3"

    local ledger_files=()
    for protein in $proteins; do
        local ledger="${task_root}/${protein}/${ledger_rel}"
        [ -f "$ledger" ] && ledger_files+=("$ledger")
    done

    if [ ${#ledger_files[@]} -eq 0 ]; then
        echo "0 0"
        return 0
    fi

    awk -F'\t' '
        NF >= 5 {
            protein = $3
            if (protein == "") {
                protein = FILENAME
            }
            key = protein ":" $4
            latest[key] = $5
        }
        END {
            success = 0
            failed = 0
            for (k in latest) {
                if (latest[k] == "success") {
                    success++
                } else if (latest[k] == "failed") {
                    failed++
                }
            }
            printf "%d %d\n", success, failed
        }
    ' "${ledger_files[@]}"
}

render_fine_method_status() {
    local task_root="$1"
    local proteins="$2"
    local ledger_rel="$3"
    local job_file_rel="$4"

    local expected success failed
    expected=$(collect_method_expected_jobs "$task_root" "$proteins" "$job_file_rel")
    read -r success failed <<< "$(collect_method_ledger_counts "$task_root" "$proteins" "$ledger_rel")"

    local processed=$((success + failed))

    if [ "$expected" -eq 0 ]; then
        echo "Pending|${success}|${failed}"
        return 0
    fi

    if [ "$processed" -eq 0 ]; then
        echo "Pending|${success}|${failed}"
        return 0
    fi

    if [ "$processed" -lt "$expected" ]; then
        local status="Processing ${processed}/${expected}"
        if [ "$failed" -gt 0 ]; then
            status="${status} (failed: ${failed})"
        fi
        echo "${status}|${success}|${failed}"
        return 0
    fi

    if [ "$failed" -gt 0 ]; then
        echo "Done (failed: ${failed})|${success}|${failed}"
    else
        echo "Done|${success}|${failed}"
    fi
}

write_progress_snapshot() {
    local task_root="${1:-$PROGRESS_MONITOR_TASK_ROOT}"
    local proteins="${2:-$PROGRESS_MONITOR_PROTEINS}"
    local pipeline_type="${3:-$PROGRESS_MONITOR_PIPELINE_TYPE}"

    if [ -z "$task_root" ] || [ -z "$proteins" ]; then
        return 0
    fi

    local status_file="${task_root}/.pipeline_progress.txt"
    local tmp_file="${status_file}.tmp"

    local total_success=0
    local total_failed=0

    local af3_line boltz2_line rf_line vina_line
    local af3_success af3_failed boltz2_success boltz2_failed rf_success rf_failed vina_success vina_failed

    IFS='|' read -r af3_line af3_success af3_failed <<< "$(render_fine_method_status \
        "$task_root" "$proteins" \
        "fine_screening/AF3/output/job_results.tsv" \
        "fine_screening/AF3/output/job_ids.txt")"
    IFS='|' read -r boltz2_line boltz2_success boltz2_failed <<< "$(render_fine_method_status \
        "$task_root" "$proteins" \
        "fine_screening/Boltz2/output/job_results.tsv" \
        "fine_screening/Boltz2/output/job_ids.txt")"
    IFS='|' read -r rf_line rf_success rf_failed <<< "$(render_fine_method_status \
        "$task_root" "$proteins" \
        "fine_screening/RoseTTAFold/protein_ligand/output/job_results.tsv" \
        "fine_screening/RoseTTAFold/protein_ligand/job_ids.txt")"
    IFS='|' read -r vina_line vina_success vina_failed <<< "$(render_fine_method_status \
        "$task_root" "$proteins" \
        "fine_screening/Vina/output/job_results.tsv" \
        "fine_screening/Vina/output/job_ids.txt")"

    total_success=$((af3_success + boltz2_success + rf_success + vina_success))
    total_failed=$((af3_failed + boltz2_failed + rf_failed + vina_failed))

    {
        echo "Status:"
        echo "--- initial screening ---"

        if [ "$pipeline_type" = "af3-first" ]; then
            echo "GraphDTA: N/A"
            echo "HMSA: N/A"
            echo "ColdDTA: N/A"
            echo "DrugLAMP: N/A"
            echo "ConPLex: N/A"
        else
            local expected_initial
            expected_initial=$(count_input_chunks "$task_root" "$proteins")
            local graphdta_done hmsa_done colddta_done druglamp_done conplex_done
            graphdta_done=$(count_initial_predictions "$task_root" "$proteins" "GraphDTA")
            hmsa_done=$(count_initial_predictions "$task_root" "$proteins" "HMSA")
            colddta_done=$(count_initial_predictions "$task_root" "$proteins" "ColdDTA")
            druglamp_done=$(count_initial_predictions "$task_root" "$proteins" "DrugLAMP")
            conplex_done=$(count_initial_predictions "$task_root" "$proteins" "ConPLex")

            echo "GraphDTA: $(format_simple_status "$graphdta_done" "$expected_initial")"
            echo "HMSA: $(format_simple_status "$hmsa_done" "$expected_initial")"
            echo "ColdDTA: $(format_simple_status "$colddta_done" "$expected_initial")"
            echo "DrugLAMP: $(format_simple_status "$druglamp_done" "$expected_initial")"
            echo "ConPLex: $(format_simple_status "$conplex_done" "$expected_initial")"
        fi

        echo ""
        echo "--- fine screening ---"
        echo "AF3: ${af3_line}"
        echo "Boltz2: ${boltz2_line}"
        echo "RoseTTAFold: ${rf_line}"
        echo "Vina: ${vina_line}"
        echo ""
        echo "Summary:"
        echo "${total_success} jobs succeeded"
        echo "${total_failed} jobs failed"
        echo ""
        echo "Updated: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
        echo "Pipeline: ${pipeline_type:-unknown}"
    } > "$tmp_file"

    mv "$tmp_file" "$status_file"
    PROGRESS_STATUS_FILE="$status_file"
}

progress_monitor_loop() {
    while true; do
        sleep "$PROGRESS_MONITOR_INTERVAL" || break
        write_progress_snapshot
    done
}

start_progress_monitor() {
    local task_root="$1"
    local proteins="$2"
    local pipeline_type="$3"
    local interval="${4:-300}"

    if [ -n "$PROGRESS_MONITOR_PID" ] && kill -0 "$PROGRESS_MONITOR_PID" 2>/dev/null; then
        stop_progress_monitor
    fi

    PROGRESS_MONITOR_TASK_ROOT="$task_root"
    PROGRESS_MONITOR_PROTEINS="$proteins"
    PROGRESS_MONITOR_PIPELINE_TYPE="$pipeline_type"
    PROGRESS_MONITOR_INTERVAL="$interval"

    write_progress_snapshot

    progress_monitor_loop &
    PROGRESS_MONITOR_PID=$!
    log_info "Started progress monitor (PID: ${PROGRESS_MONITOR_PID})"
}

refresh_progress_monitor() {
    write_progress_snapshot
}

stop_progress_monitor() {
    if [ -n "$PROGRESS_MONITOR_PID" ] && kill -0 "$PROGRESS_MONITOR_PID" 2>/dev/null; then
        kill "$PROGRESS_MONITOR_PID" 2>/dev/null || true
        wait "$PROGRESS_MONITOR_PID" 2>/dev/null || true
    fi

    if [ -n "$PROGRESS_MONITOR_TASK_ROOT" ]; then
        write_progress_snapshot
    fi

    PROGRESS_MONITOR_PID=""
}
