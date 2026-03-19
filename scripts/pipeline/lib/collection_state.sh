#!/bin/bash
#
# Collection State Management Library
# Manages JSON state files for tracking incremental results collection
#

# Initialize collection state file
# Args: protein, method, expected_count, state_file
init_collection_state() {
    local protein="$1"
    local method="$2"
    local expected_count="$3"
    local state_file="$4"

    # Create directory if needed
    mkdir -p "$(dirname "$state_file")"

    cat > "$state_file" <<EOF
{
  "method": "${method}",
  "protein": "${protein}",
  "last_collection": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "total_expected": ${expected_count},
  "collected_compounds": {},
  "summary_csv_rows": 0,
  "last_warning": null
}
EOF
}

# Mark a compound as collected
# Args: state_file, compound_id, has_data, [warning_msg]
mark_compound_collected() {
    local state_file="$1"
    local compound_id="$2"
    local has_data="$3"
    local warning_msg="${4:-}"

    if [ ! -f "$state_file" ]; then
        echo "Error: State file not found: $state_file" >&2
        return 1
    fi

    local timestamp
    timestamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)

    # Build JSON entry
    local entry="{\"timestamp\": \"${timestamp}\", \"status\": \"collected\", \"has_data\": ${has_data}"
    if [ -n "$warning_msg" ]; then
        entry="${entry}, \"warning\": \"${warning_msg}\""
    fi
    entry="${entry}}"

    # Update state file using jq
    jq --arg id "$compound_id" --argjson entry "$entry" \
        '.collected_compounds[$id] = $entry | .last_collection = "'$timestamp'"' \
        "$state_file" > "${state_file}.tmp" && mv "${state_file}.tmp" "$state_file"
}

# Check if a compound has been collected
# Args: state_file, compound_id
# Returns: 0 if collected, 1 if not
is_compound_collected() {
    local state_file="$1"
    local compound_id="$2"

    if [ ! -f "$state_file" ]; then
        return 1
    fi

    jq -e --arg id "$compound_id" '.collected_compounds | has($id)' "$state_file" > /dev/null 2>&1
}

# Get list of uncollected compound IDs
# Args: state_file, all_compounds_file (one ID per line)
# Outputs: List of uncollected IDs (one per line)
get_uncollected_compounds() {
    local state_file="$1"
    local all_compounds_file="$2"

    if [ ! -f "$state_file" ]; then
        # No state file = all uncollected
        cat "$all_compounds_file"
        return 0
    fi

    # Get collected IDs
    local collected_ids
    collected_ids=$(jq -r '.collected_compounds | keys[]' "$state_file")

    # Filter out collected from all
    while IFS= read -r compound_id; do
        if ! echo "$collected_ids" | grep -qx "$compound_id"; then
            echo "$compound_id"
        fi
    done < "$all_compounds_file"
}

# Get collection statistics
# Args: state_file
# Outputs: JSON with stats
get_collection_stats() {
    local state_file="$1"

    if [ ! -f "$state_file" ]; then
        echo '{"collected": 0, "expected": 0, "percentage": 0}'
        return 0
    fi

    jq '{
        collected: (.collected_compounds | length),
        expected: .total_expected,
        percentage: ((.collected_compounds | length) * 100 / .total_expected),
        last_collection: .last_collection,
        has_warnings: ([.collected_compounds[] | select(.has_data == false)] | length > 0)
    }' "$state_file"
}

# Update collection timestamp
# Args: state_file
update_collection_timestamp() {
    local state_file="$1"

    if [ ! -f "$state_file" ]; then
        return 1
    fi

    local timestamp
    timestamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)

    jq --arg ts "$timestamp" '.last_collection = $ts' "$state_file" > "${state_file}.tmp" \
        && mv "${state_file}.tmp" "$state_file"
}

# Get compounds with warnings (missing data)
# Args: state_file
# Outputs: List of compound IDs with warnings
get_compounds_with_warnings() {
    local state_file="$1"

    if [ ! -f "$state_file" ]; then
        return 0
    fi

    jq -r '.collected_compounds | to_entries[] | select(.value.has_data == false) | .key' "$state_file"
}

# Update summary CSV row count
# Args: state_file, row_count
update_summary_row_count() {
    local state_file="$1"
    local row_count="$2"

    if [ ! -f "$state_file" ]; then
        return 1
    fi

    jq --arg count "$row_count" '.summary_csv_rows = ($count | tonumber)' "$state_file" > "${state_file}.tmp" \
        && mv "${state_file}.tmp" "$state_file"
}
