#!/bin/bash
#
# Configuration Loader Library
# Load pipeline configuration from YAML files
#

# Source dependencies (only if not already sourced)
if [ -z "$(type -t log_info)" ]; then
    LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${LIB_DIR}/logger.sh"
fi

# Parse YAML to JSON using Python
parse_yaml_to_json() {
    local yaml_file="$1"

    if [ ! -f "$yaml_file" ]; then
        log_error "YAML file not found: ${yaml_file}"
        return 1
    fi

    python3 -c "
import yaml
import json
import sys

try:
    with open('${yaml_file}', 'r') as f:
        data = yaml.safe_load(f)
    print(json.dumps(data))
except Exception as e:
    print(f'Error parsing YAML: {e}', file=sys.stderr)
    sys.exit(1)
" 2>/dev/null

    if [ $? -ne 0 ]; then
        log_error "Failed to parse YAML file: ${yaml_file}"
        return 1
    fi
}

# Load configuration from YAML file
load_yaml_config() {
    local yaml_file="$1"

    log_info "Loading configuration from: ${yaml_file}"

    # Check if PyYAML is available
    if ! python3 -c "import yaml" 2>/dev/null; then
        log_error "PyYAML not installed. Install with: pip install pyyaml"
        return 1
    fi

    # Parse YAML to JSON
    local config_json
    config_json=$(parse_yaml_to_json "$yaml_file")

    if [ $? -ne 0 ]; then
        return 1
    fi

    echo "$config_json"
}

# Extract value from config JSON
get_config_value() {
    local config_json="$1"
    local key_path="$2"
    local default_value="${3:-}"

    local value
    value=$(echo "$config_json" | jq -r "$key_path" 2>/dev/null)

    if [ "$value" = "null" ] || [ -z "$value" ]; then
        echo "$default_value"
    else
        echo "$value"
    fi
}

# Extract array from config JSON
get_config_array() {
    local config_json="$1"
    local key_path="$2"

    echo "$config_json" | jq -r "${key_path}[]" 2>/dev/null
}

# Check if key exists in config
config_has_key() {
    local config_json="$1"
    local key_path="$2"

    local value
    value=$(echo "$config_json" | jq -r "$key_path" 2>/dev/null)

    if [ "$value" = "null" ]; then
        return 1
    else
        return 0
    fi
}
