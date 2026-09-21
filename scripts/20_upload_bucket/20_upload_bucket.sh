#!/usr/bin/env bash
# Stage 2 — Upload the Stage 1 freeze to the staging bucket, step by step.
#
# Each numbered subfolder (NN_*/) under scripts/20_upload_bucket/ exposes a single
# build.sh entry point. This driver runs them in numeric order; the FOLDER SEQUENCE below
# is the run order. The shape is Stage 1's — same STEP_DIRS list, same STEPS filter, same
# per-step log and report TSV — because a stage that reads like the other stages is one
# fewer thing to learn.
#
# A step's build.sh exit code: 0 = PASS, 77 = SKIP, anything else = FAIL.
#
# Promotion to production is intentionally NOT here — see `make promote`.
#
# Env knobs:
#   STEPS="01 03"    run only steps whose folder name contains one of these tokens
#   STOP_ON_FAIL=... abort on the first FAIL (default: true)
#   APPLY=1          step 04 actually uploads (default: dry run, recorded as SKIP)
#   SEED_STAGING=1   step 02 mirrors production into staging (once per release cycle)
#   FREEZE_DIR=...   read the freeze from somewhere other than ${PRECOVID_ROOT}/staging/freeze
#   EXECUTOR=slurm   submit each step as its own SLURM job (config/slurm.json). The
#                    compute node needs the same gcloud credentials as the submit
#                    host — on SCG ~/.config/gcloud is on shared storage, so it does.
#
# Output: logs/upload_<step>.log per step; logs/upload_bucket_report.tsv summary.

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

STEPS="${STEPS:-}"
STOP_ON_FAIL="${STOP_ON_FAIL:-true}"
readonly SKIP_RC=77

STAGE_DIR="${PIPELINE_ROOT}/scripts/20_upload_bucket"
REPORT_TSV="${LOG_DIR}/upload_bucket_report.tsv"
printf 'status\tstep\tdetail\n' > "${REPORT_TSV}"

log "Stage 2 — upload to bucket (STOP_ON_FAIL=${STOP_ON_FAIL})"
log "staging bucket: ${STAGING_BUCKET}"
[[ -n "${STEPS}" ]] && log "Step filter: ${STEPS}"
if [[ "${APPLY:-0}" == "1" ]]; then
  warn "APPLY=1 — step 04 will WRITE to ${STAGING_BUCKET}"
else
  log "dry run — step 04 will plan only (APPLY=1 to upload)"
fi

# Gate: tooling, bucket read access, the freeze, and the manifest's input object.
bash "${STAGE_DIR}/check_required_inputs.sh" \
  || die "Stage 2 aborted: required inputs missing (see above)."

# The ordered run sequence. Adding/renaming a step = editing this list + the folder.
STEP_DIRS=(
  01_snapshot_buckets
  02_seed_staging
  03_diff_freeze_vs_staging
  04_upload_to_staging
  05_validate_structure
)

# step_selected DIR -> 0 if it should run given the STEPS filter
step_selected() {
  [[ -z "${STEPS}" ]] && return 0
  local want; for want in ${STEPS}; do [[ "$1" == *"${want}"* ]] && return 0; done
  return 1
}

ran=0; passed=0; failed=0; skipped=0
for step in "${STEP_DIRS[@]}"; do
  step_selected "${step}" || continue
  bs="${STAGE_DIR}/${step}/build.sh"
  if [[ ! -f "${bs}" ]]; then
    printf 'FAIL\t%s\tbuild.sh not found\n' "${step}" >> "${REPORT_TSV}"
    err "${step}: build.sh not found"; failed=$((failed+1))
    [[ "${STOP_ON_FAIL}" == "true" ]] && die "stopping (STOP_ON_FAIL) at ${step}"
    continue
  fi

  log "step ${step}"
  logfile="${LOG_DIR}/upload_${step}.log"
  ran=$((ran+1))
  set +e
  run_step "20_upload_bucket/${step}" "${logfile}" bash "${bs}"
  rc=$?
  set -e

  tail_msg="$(log_tail_detail "${logfile}")"
  if (( rc == 0 )); then
    printf 'PASS\t%s\tok\n' "${step}" >> "${REPORT_TSV}"
    ok "${step} — see $(basename "${logfile}")"; passed=$((passed+1))
  elif (( rc == SKIP_RC )); then
    printf 'SKIP\t%s\t%s\n' "${step}" "${tail_msg}" >> "${REPORT_TSV}"
    warn "${step} — SKIP (${tail_msg})"; skipped=$((skipped+1))
  else
    printf 'FAIL\t%s\trc=%d: %s\n' "${step}" "${rc}" "${tail_msg}" >> "${REPORT_TSV}"
    err "${step} failed (rc=${rc}) — see ${logfile}"; failed=$((failed+1))
    [[ "${STOP_ON_FAIL}" == "true" ]] && die "stopping (STOP_ON_FAIL) at ${step}"
  fi
done

log "Stage 2: ${ran} ran, ${passed} PASS, ${skipped} SKIP, ${failed} FAIL — see ${REPORT_TSV}"
(( failed == 0 )) || die "Stage 2 finished with ${failed} failed step(s)."
# Only step 05 writes structure_validation_*.tsv, and it exits SKIP while step 03 still
# lists files waiting to be uploaded — the default posture, since step 04 needs APPLY=1.
# Naming that report unconditionally would send a dry run to read the previous release's
# file and take this freeze for validated.
if (( skipped > 0 )); then
  ok "Stage 2 complete — ${skipped} step(s) skipped, nothing written to GCS. The upload plan is in logs/upload_plan_*.tsv; re-run with APPLY=1 to upload and validate."
else
  ok "Stage 2 complete. Review logs/structure_validation_*.tsv before promote."
fi
