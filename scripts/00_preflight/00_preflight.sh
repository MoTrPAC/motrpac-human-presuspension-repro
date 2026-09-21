#!/usr/bin/env bash
# Stage 0 — Preflight: verify permissions, tooling, repos and external assets
# BEFORE any expensive regeneration runs. Fails fast on missing REQUIRED prereqs;
# WARNs (does not fail) on optional assets. Consortium access is required: every
# node runs, so an unreadable gated bucket is a FAIL, not a degraded run.
#
# Output: logs/preflight_report.tsv (status<TAB>check<TAB>detail) + console summary.
# Exit:   non-zero if any REQUIRED check FAILs.

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

REPORT_TSV="${LOG_DIR}/preflight_report.tsv"
report_init "${REPORT_TSV}"
log "Preflight — report: ${REPORT_TSV}"

# ---- 1. Required tools on PATH ---------------------------------------------
for tool in "${GSUTIL}" "${GCLOUD}" "${RSCRIPT}" curl; do
  if have "${tool}"; then record_check PASS "tool:${tool}" "found at $(command -v "${tool}")"
  else record_check FAIL "tool:${tool}" "not on PATH"; fi
done

# ---- 2. Required R packages ------------------------------------------------
# Parse Imports from both DESCRIPTIONs would be ideal; here we check the load-
# bearing Bioconductor/CRAN deps that block the build if absent.
if have "${RSCRIPT}"; then
  R_PKGS=(dplyr magrittr data.table tibble tidyr ggplot2 jsonlite here devtools \
          Biobase ComplexHeatmap Mfuzz TMSig variancePartition MotrpacBicQC)
  # Step 16 only. run_SCION() used to check these itself with check_package_installation()
  # at call time; the file that did is gone, and discovering a missing randomForest after
  # the matrices have loaded costs the whole step. e1071 supplies cmeans(), which
  # Mfuzz::mfuzz() calls off the search path rather than through its own namespace.
  if [[ "$(printf '%s' "${RUN_SCION:-TRUE}" | tr '[:lower:]' '[:upper:]')" == "TRUE" ]]; then
    R_PKGS+=(parallel doParallel randomForest e1071)
  fi
  missing_pkgs="$("${RSCRIPT}" -e '
    pkgs <- commandArgs(trailingOnly=TRUE)
    miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
    cat(paste(miss, collapse=" "))
  ' "${R_PKGS[@]}" 2>/dev/null || echo "R_INVOCATION_FAILED")"
  if [[ "${missing_pkgs}" == "R_INVOCATION_FAILED" ]]; then
    record_check FAIL "r_packages" "could not invoke Rscript to check packages"
  elif [[ -z "${missing_pkgs}" ]]; then
    record_check PASS "r_packages" "all ${#R_PKGS[@]} required packages installed"
  else
    record_check FAIL "r_packages" "missing: ${missing_pkgs}"
  fi
fi

# ---- 3. gcloud authentication ----------------------------------------------
if have "${GCLOUD}"; then
  # ADC is only a hard requirement if bucket access fails (§5 is the real gate);
  # gsutil itself authenticates off the gcloud login, so WARN not FAIL here.
  if "${GCLOUD}" auth application-default print-access-token >/dev/null 2>&1; then
    record_check PASS "gcloud_adc" "application-default credentials active"
  else
    record_check WARN "gcloud_adc" "no ADC — fine if bucket checks pass; else run: gcloud auth application-default login"
  fi
  active_acct="$("${GCLOUD}" auth list --filter=status:ACTIVE --format='value(account)' 2>/dev/null | head -n1)"
  [[ -n "${active_acct}" ]] \
    && record_check PASS "gcloud_account" "active: ${active_acct}" \
    || record_check WARN "gcloud_account" "no active gcloud account"
fi

# ---- 4. In-scope repos present ---------------------------------------------
for repo in "${DATA_PKG_REPO}" "${ANALYSIS_PKG_REPO}"; do
  [[ -d "${repo}" ]] \
    && record_check PASS "repo:$(basename "${repo}")" "present" \
    || record_check FAIL "repo:$(basename "${repo}")" "not found at ${repo}"
done

# ---- 5. Bucket access (consortium read is required) ------------------------
if have "${GSUTIL}"; then
  if gcs_can_read "${PRODUCTION_BUCKET}"; then
    record_check PASS "bucket_read" "readable: ${PRODUCTION_BUCKET}"
  else
    record_check FAIL "bucket_read" "cannot read ${PRODUCTION_BUCKET} — consortium access is required; request it or run: gcloud auth login"
  fi
  if gcs_can_write "${STAGING_BUCKET}"; then
    record_check PASS "bucket_write" "writable: ${STAGING_BUCKET}"
  else
    record_check WARN "bucket_write" "cannot write ${STAGING_BUCKET} (Stage 2 upload will be blocked)"
  fi
fi

# ---- 6. External downloadable assets (from external_assets.tsv) -------------
# Resolve a manifest `file` location and echo the path that exists: absolute paths
# as-is, then relative to PIPELINE_ROOT (every non-gated raw source is vendored in
# this repo), then ~/Downloads for the assets a provider only hands out through a
# logged-in browser. Non-zero if none of them exist.
resolve_asset_file() {
  local loc="$1" cand
  for cand in "${loc}" "${PIPELINE_ROOT}/${loc}" "${HOME}/Downloads/$(basename "${loc}")"; do
    if [[ -e "${cand}" ]]; then printf '%s' "${cand}"; return 0; fi
  done
  return 1
}

ASSETS_TSV="${PIPELINE_ROOT}/scripts/00_preflight/config/external_assets.tsv"
if [[ -f "${ASSETS_TSV}" ]]; then
  # columns: id kind location consumed_by gated notes  (tab-separated, 1 header row)
  while IFS=$'\t' read -r id kind location consumed_by gated notes; do
    [[ "${id}" == "id" || -z "${id}" ]] && continue
    case "${kind}" in
      url)
        if url_reachable "${location}"; then record_check PASS "asset:${id}" "reachable"
        else record_check WARN "asset:${id}" "unreachable: ${location}"; fi ;;
      file)
        if [[ "${location}" == gs://* ]]; then
          # Single gated object: stat it now rather than deferring to Stage 2 —
          # §5 only proves the enclosing bucket lists, not that this file is there.
          if ! have "${GSUTIL}"; then
            record_check WARN "asset:${id}" "cannot stat without gsutil: ${location}"
          elif gcs_can_read "${location}"; then
            record_check PASS "asset:${id}" "readable on bucket"
          else
            record_check WARN "asset:${id}" "not readable on bucket: ${location}"
          fi
        elif found="$(resolve_asset_file "${location}")"; then
          record_check PASS "asset:${id}" "present: ${found#"${PIPELINE_ROOT}"/}"
        else
          record_check WARN "asset:${id}" "local file not found: ${location}"
        fi ;;
      bucket|auth) : ;;  # covered by sections 3 and 5
      *) record_check WARN "asset:${id}" "unknown kind '${kind}'" ;;
    esac
  done < "${ASSETS_TSV}"
else
  record_check WARN "assets_manifest" "external_assets.tsv not found"
fi

# ---- 7. GMT files for MOLECULAR_SIGNATURES ---------------------------------
# Check the 7 GMTs Stage 1 step 03 consumes, one PASS/WARN row each.
#
# The canonical read dir is THIS stage's output/: the vendored
# 03_build_molecular_signatures/MOLECULAR_SIGNATURES.R reads it instead of the
# Analysis package's data-raw/gmt_files/, so that package is not consulted here.
# data-raw/build_data_objects.sh populates output/ — building five GMTs from the
# vendored raw sources and copying the two download-only MSigDB ones verbatim.
#
# An unbuilt GMT is a WARN, not a FAIL, because `make data-objects` produces it.
# The detail line distinguishes the two cases that need different fixes: source
# vendored (just build) vs source absent (re-download it first).
#   read_dir  : where the vendored MOLECULAR_SIGNATURES.R reads (output/)
#   vendor_dir: raw inputs — sources/ per-builder, prebuilt/ for download-only
check_gmt_files() {
  local stage="${PIPELINE_ROOT}/scripts/00_preflight"
  local read_dir="${stage}/output"
  local vendor_dir="${stage}/data-raw/gmt_processing"
  local tsv="${stage}/config/gmt_files.tsv"
  [[ -f "${tsv}" ]] || { record_check WARN "gmt_manifest" "${tsv} not found"; return 0; }

  local unbuilt=0 src
  while IFS=$'\t' read -r gmt_file source_file url builder license; do
    [[ "${gmt_file}" == "gmt_file" || -z "${gmt_file}" ]] && continue
    if [[ -f "${read_dir}/${gmt_file}" ]]; then
      record_check PASS "gmt:${gmt_file}" "built — output/${gmt_file}"
      continue
    fi
    unbuilt=1
    # A parenthesized source_file ("(downloaded pre-built)", "(RefMet API query)")
    # means there is no local raw input to process: the artifact itself is vendored
    # under prebuilt/, either download-only or a frozen API snapshot.
    if [[ "${source_file}" == \(* ]]; then
      src="${vendor_dir}/prebuilt/${gmt_file}"
    else
      src="${vendor_dir}/sources/${source_file}"
    fi
    if [[ -f "${src}" ]]; then
      record_check WARN "gmt:${gmt_file}" "not built; source vendored (${src#"${PIPELINE_ROOT}"/}) — run: make data-objects"
    else
      record_check WARN "gmt:${gmt_file}" "not built and source missing (${src#"${PIPELINE_ROOT}"/}) — build via ${builder}; download source from ${url} (${license})"
    fi
  done < "${tsv}"

  if (( unbuilt )); then
    warn "One or more GMTs absent from scripts/00_preflight/output/ — run: make data-objects"
  fi
  return 0
}

check_gmt_files

# ---- 8. SLURM submission (only when the run intends to use it) --------------
# Silent on a laptop: with EXECUTOR unset and config/slurm.json disabled there is
# nothing to verify. When a run WILL submit, the account and the submit tools are
# checked here rather than at the first sbatch, because discovering an unset account
# 40 minutes into stage 1 costs the whole stage.
if slurm_enabled; then
  for tool in sbatch squeue scancel; do
    if have "${tool}"; then record_check PASS "slurm:${tool}" "found at $(command -v "${tool}")"
    else record_check FAIL "slurm:${tool}" "EXECUTOR=slurm but ${tool} is not on PATH"; fi
  done

  slurm_account="$(slurm_cfg data account)"
  if [[ -z "${slurm_account}" || "${slurm_account}" == "CHANGE_ME" ]]; then
    record_check FAIL "slurm_account" "set defaults.account in config/slurm.json (see docs/SLURM.md)"
  elif have sacctmgr && ! sacctmgr -nP show associations user="${USER:-$(id -un)}" account="${slurm_account}" 2>/dev/null | grep -q .; then
    # An account you are not associated with is accepted by sbatch and then rejected
    # by the scheduler, so the job never starts and nothing says why.
    record_check WARN "slurm_account" "${USER:-$(id -un)} has no association with account '${slurm_account}'"
  else
    record_check PASS "slurm_account" "${slurm_account}"
  fi

  slurm_partition="$(slurm_cfg data partition)"
  if [[ -n "${slurm_partition}" ]] && have sinfo; then
    sinfo -h -p "${slurm_partition}" -o '%P' 2>/dev/null | grep -q . \
      && record_check PASS "slurm_partition" "${slurm_partition}" \
      || record_check FAIL "slurm_partition" "partition '${slurm_partition}' unknown to this cluster"
  fi

  record_check PASS "slurm_config" "$(basename "${SLURM_CONFIG}") loaded — stage 1 asks for cpus=$(slurm_cfg data cpus_per_task), mem=$(slurm_cfg data mem), time=$(slurm_cfg data time)"
fi

# ---- 9. Summarize ---------------------------------------------------------
echo
log "Preflight summary: ${CHECK_FAILS} FAIL, ${CHECK_WARNS} WARN."
if (( CHECK_FAILS > 0 )); then
  die "Preflight failed on ${CHECK_FAILS} required check(s). See ${REPORT_TSV}."
fi
ok "Preflight passed."
