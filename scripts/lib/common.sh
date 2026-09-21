#!/usr/bin/env bash
# Shared helpers for every stage script. Source this first:
#   source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
# It sets PIPELINE_ROOT, loads config/pipeline.env, and defines logging +
# gsutil/R helpers used across stages.

set -euo pipefail

# ---- Locate the project root and load config -------------------------------
# lib/ lives at scripts/lib/, so the project root is two levels up.
_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${_LIB_DIR}/../.." && pwd)"
export PIPELINE_ROOT

# shellcheck source=/dev/null
source "${PIPELINE_ROOT}/config/pipeline.env"

LOG_DIR="${PIPELINE_ROOT}/logs"
mkdir -p "${LOG_DIR}"

# ---- Logging ---------------------------------------------------------------
_c_reset=$'\033[0m'; _c_red=$'\033[31m'; _c_grn=$'\033[32m'; _c_yel=$'\033[33m'; _c_blu=$'\033[34m'
log()  { printf '%s[%s]%s %s\n' "${_c_blu}" "$(_ts)" "${_c_reset}" "$*"; }
ok()   { printf '%s[ OK ]%s %s\n' "${_c_grn}" "${_c_reset}" "$*"; }
warn() { printf '%s[WARN]%s %s\n' "${_c_yel}" "${_c_reset}" "$*" >&2; }
err()  { printf '%s[FAIL]%s %s\n' "${_c_red}" "${_c_reset}" "$*" >&2; }
die()  { err "$*"; exit 1; }

# Timestamp helper (no subshell-in-format surprises).
_ts() { date +'%Y-%m-%d %H:%M:%S'; }

# ---- Check-result accounting -----------------------------------------------
# Stage scripts append rows to a report TSV via record_check.
CHECK_FAILS=0
CHECK_WARNS=0
report_init() {
  REPORT_TSV="$1"
  printf 'status\tcheck\tdetail\n' > "${REPORT_TSV}"
}
# record_check STATUS CHECK DETAIL   (STATUS in PASS|WARN|FAIL)
record_check() {
  local status="$1" check="$2" detail="${3:-}"
  printf '%s\t%s\t%s\n' "${status}" "${check}" "${detail}" >> "${REPORT_TSV}"
  case "${status}" in
    PASS) ok   "${check} — ${detail}" ;;
    WARN) warn "${check} — ${detail}"; CHECK_WARNS=$((CHECK_WARNS+1)) ;;
    FAIL) err  "${check} — ${detail}"; CHECK_FAILS=$((CHECK_FAILS+1)) ;;
  esac
}

# ---- Small utilities -------------------------------------------------------
have() { command -v "$1" >/dev/null 2>&1; }

# gsutil read/write probes (return 0 on success, non-fatal).
gcs_can_read()  { "${GSUTIL}" ls "$1" >/dev/null 2>&1; }
gcs_can_write() {
  local bucket="$1" probe="$1/.precovid_repro_write_probe_$$"
  if printf 'probe' | "${GSUTIL}" cp - "${probe}" >/dev/null 2>&1; then
    "${GSUTIL}" rm "${probe}" >/dev/null 2>&1 || true
    return 0
  fi
  return 1
}

# URL reachability (HEAD; falls back to a range GET for servers that reject HEAD).
url_reachable() {
  local url="$1"
  curl -fsS -m 15 -I "${url}" >/dev/null 2>&1 && return 0
  curl -fsS -m 15 -r 0-0 "${url}" >/dev/null 2>&1
}

# Run an R script inside a given repo, tee-ing output to a log.
run_r() {
  local repo="$1" rscript="$2" logfile="$3"
  ( cd "${repo}" && "${RSCRIPT}" "${rscript}" ) 2>&1 | tee "${logfile}"
}

# The detail column of a step's report row: the last line the step actually
# printed. Strips ANSI colour, tabs (which would split the TSV), and the "[slurm]"
# markers the batch wrapper adds — those are about the job, not about the step, and
# a failed step's own last line is what makes the report readable.
log_tail_detail() {
  grep -v '^\[slurm\]' "$1" 2>/dev/null \
    | tail -n1 \
    | sed $'s/\033\\[[0-9;]*m//g' \
    | tr '\t' ' '
}

# Map a manifest repo key ("analysis" | "data") to its filesystem path.
resolve_repo() {
  case "$1" in
    analysis) printf '%s' "${ANALYSIS_PKG_REPO}" ;;
    data)     printf '%s' "${DATA_PKG_REPO}" ;;
    *) return 1 ;;
  esac
}

# ---- Execution backend ------------------------------------------------------
# Defines run_step, which the stage drivers use instead of calling build.sh
# directly. Local by default; submits to SLURM when EXECUTOR=slurm or
# config/slurm.json sets enabled=true. See docs/SLURM.md.
# shellcheck source=/dev/null
source "${_LIB_DIR}/slurm.sh"
