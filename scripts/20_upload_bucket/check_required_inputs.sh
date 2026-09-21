#!/usr/bin/env bash
# Stage 2 required-inputs gate — run at the START of 20_upload_bucket, before any step.
# Same posture as the Stage 1 gate: fail fast HERE with one consolidated list rather than
# have a step die part-way through against a bucket it cannot read. HARD gate.
#
# Scope = what must PRE-EXIST for the stage to run, not what the stage produces:
#   - the Stage 1 freeze (this stage uploads it)
#   - gcloud/gsutil on PATH, authenticated, with read access to both buckets
#   - the preflight OME_TISSUE_CODE object the manifest is generated from
#   - base64enc, for converting gsutil's base64 md5 to the hex tools::md5sum returns
#
# Write access to the staging bucket is NOT checked here. Only step 04 writes, only under
# APPLY=1, and it probes for itself — a dry run and a validation pass are useful without it.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

report_init "${LOG_DIR}/upload_required_inputs_report.tsv"
log "Stage 2 required-inputs check"

# --- Tooling ----------------------------------------------------------------
for tool in "${GSUTIL}" "${GCLOUD}" "${RSCRIPT}"; do
  if have "${tool}"; then
    record_check PASS "tool:${tool}" "on PATH"
  else
    record_check FAIL "tool:${tool}" "not found on PATH"
  fi
done

if "${RSCRIPT}" -e 'quit(status = !requireNamespace("base64enc", quietly = TRUE))' 2>/dev/null; then
  record_check PASS "rpkg:base64enc" "installed"
else
  record_check FAIL "rpkg:base64enc" "not installed (install.packages(\"base64enc\"))"
fi

# --- Bucket access ----------------------------------------------------------
# Read-only probes. The buckets themselves are pinned in config/pipeline.env.
for entry in "production::${PRODUCTION_BUCKET}" "staging::${STAGING_BUCKET}"; do
  label="${entry%%::*}"; bucket="${entry##*::}"
  if gcs_can_read "${bucket}"; then
    record_check PASS "bucket:${label}" "readable at ${bucket}"
  else
    record_check FAIL "bucket:${label}" "not readable: ${bucket} (re-auth: gcloud auth login)"
  fi
done

# --- Manifest input ---------------------------------------------------------
ome_tissue="${PIPELINE_ROOT}/scripts/00_preflight/data/OME_TISSUE_CODE.rds"
if [[ -f "${ome_tissue}" ]]; then
  record_check PASS "preflight:OME_TISSUE_CODE" "present"
else
  record_check FAIL "preflight:OME_TISSUE_CODE" "MISSING: ${ome_tissue#"${PIPELINE_ROOT}"/} (run Stage 0: make preflight data-objects)"
fi

# --- The freeze -------------------------------------------------------------
# What this stage uploads. FREEZE_DIR overrides the location, which is how a checkout
# without its own staging/ (it is gitignored, so a worktree has none) runs against a
# freeze built elsewhere.
freeze_dir="${FREEZE_DIR:-${STAGING_DIR}/freeze}"
if [[ -d "${freeze_dir}" ]]; then
  n="$(find "${freeze_dir}" -type f \( -name '*.txt' -o -name '*.txt.gz' -o -name '*.csv' -o -name '*.html' \) | wc -l | tr -d ' ')"
  if (( n > 0 )); then
    record_check PASS "freeze" "${n} upload artifact(s) under ${freeze_dir}"
  else
    record_check FAIL "freeze" "no upload artifacts under ${freeze_dir} (run Stage 1: make data)"
  fi
else
  record_check FAIL "freeze" "MISSING: ${freeze_dir} (run Stage 1: make data, or set FREEZE_DIR)"
fi

if (( CHECK_FAILS > 0 )); then
  die "${CHECK_FAILS} required input(s) missing — see ${REPORT_TSV}. Fix the above before running Stage 2."
fi
ok "All required Stage 2 inputs present (${CHECK_WARNS} warning(s))."
