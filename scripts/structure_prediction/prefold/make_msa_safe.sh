#!/bin/bash

set -euo pipefail

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*" >&2
}

log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN: $*" >&2
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

if [ "$#" -lt 5 ]; then
    log_error "Usage: $0 <in_fasta> <out_dir> <cpu> <mem> <db_templ> [extra args]"
    exit 2
fi

IN_FASTA="$1"
OUT_DIR="$2"
RFAA_ROOT="${RFAA_ROOT:-$(pwd)}"
UPSTREAM_MAKE_MSA="${RFAA_ROOT}/make_msa.sh"
MAX_ATTEMPTS=2

if [ ! -x "$UPSTREAM_MAKE_MSA" ]; then
    log_error "Upstream make_msa.sh not found or not executable: ${UPSTREAM_MAKE_MSA}"
    exit 1
fi

mkdir -p "$OUT_DIR"

cleanup_partial_outputs() {
    rm -f \
        "${OUT_DIR}/t000_.ss2" \
        "${OUT_DIR}/t000_.hhr" \
        "${OUT_DIR}/t000_.atab" \
        "${OUT_DIR}/t000_.msa0.ss2.a3m"
}

validate_ss2() {
    local ss2_file=$1
    local ss_pred
    local ss_conf

    [ -s "$ss2_file" ] || return 1
    grep -q '^>ss_pred$' "$ss2_file" || return 1
    grep -q '^>ss_conf$' "$ss2_file" || return 1

    ss_pred="$(awk '/^>ss_pred$/ {getline; print; exit}' "$ss2_file" | tr -d '[:space:]')"
    ss_conf="$(awk '/^>ss_conf$/ {getline; print; exit}' "$ss2_file" | tr -d '[:space:]')"

    [ -n "$ss_pred" ] && [ -n "$ss_conf" ]
}

validate_outputs() {
    local ok=0

    if [ ! -s "${OUT_DIR}/t000_.msa0.a3m" ]; then
        log_error "Missing or empty output: ${OUT_DIR}/t000_.msa0.a3m"
        ok=1
    fi

    if ! validate_ss2 "${OUT_DIR}/t000_.ss2"; then
        log_error "Malformed or empty output: ${OUT_DIR}/t000_.ss2"
        ok=1
    fi

    if [ ! -s "${OUT_DIR}/t000_.hhr" ]; then
        log_error "Missing or empty output: ${OUT_DIR}/t000_.hhr"
        ok=1
    fi

    if [ ! -s "${OUT_DIR}/t000_.atab" ]; then
        log_error "Missing or empty output: ${OUT_DIR}/t000_.atab"
        ok=1
    fi

    if [ "$ok" -ne 0 ]; then
        return 1
    fi
    return 0
}

run_upstream_once() {
    local attempt=$1
    shift
    local tmp_parent="${OUT_DIR}/.msa_tmp"
    local tmp_work_dir
    local rc

    mkdir -p "$tmp_parent"
    tmp_work_dir="$(mktemp -d "${tmp_parent}/run_${SLURM_JOB_ID:-nojid}_${attempt}_XXXXXX")"
    rc=0

    log_info "Attempt ${attempt}/${MAX_ATTEMPTS}: running make_msa.sh in isolated workdir ${tmp_work_dir}"

    (
        cd "$tmp_work_dir"
        "$UPSTREAM_MAKE_MSA" "$@"
    ) || rc=$?

    rm -rf "$tmp_work_dir"
    return "$rc"
}

attempt=1
while [ "$attempt" -le "$MAX_ATTEMPTS" ]; do
    if run_upstream_once "$attempt" "$@"; then
        if validate_outputs; then
            log_info "Validated MSA/template artifacts for ${IN_FASTA} in ${OUT_DIR}"
            exit 0
        fi
        log_error "Artifact validation failed for ${IN_FASTA} (attempt ${attempt}/${MAX_ATTEMPTS})"
    else
        log_error "Upstream make_msa.sh failed for ${IN_FASTA} (attempt ${attempt}/${MAX_ATTEMPTS})"
    fi

    if [ "$attempt" -lt "$MAX_ATTEMPTS" ]; then
        log_warn "Cleaning partial outputs before retry"
        cleanup_partial_outputs
    fi

    attempt=$((attempt + 1))
done

log_error "MSA/template generation failed after ${MAX_ATTEMPTS} attempts for ${IN_FASTA}"
if [ -f "${OUT_DIR}/log/make_ss.stderr" ]; then
    log_error "Inspect PSIPRED log: ${OUT_DIR}/log/make_ss.stderr"
fi
exit 1
