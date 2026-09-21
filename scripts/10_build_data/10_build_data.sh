#!/usr/bin/env bash
# Stage 1 — Build data objects, step by step.
#
# Each numbered subfolder (NN_*/) under scripts/10_build_data/ exposes a single
# build.sh entry point. This driver runs them in numeric order; the FOLDER SEQUENCE
# below is the build order (dependency/creation order). Each build.sh does one step:
# a staging copy, an adapted local R build, the qc-norm stem loop, or a not-yet-
# adapted placeholder. docs/data_objects.tsv is now the object inventory /
# regeneration ledger — it no longer drives the build.
#
# A step's build.sh exit code: 0 = PASS, 77 = SKIP (placeholder / not yet adapted),
# anything else = FAIL.
#
# Env knobs:
#   STEPS="03 06"    run only steps whose folder name contains one of these tokens
#                    (default: all)
#   STOP_ON_FAIL=... abort on the first FAIL (default: true)
#   EXECUTOR=slurm   submit each step as its own SLURM job (resources per step from
#                    config/slurm.json) and block until it finishes, instead of
#                    running it here. The report below is identical either way.
#
# Output: logs/data_<step>.log per step; logs/build_data_report.tsv summary.

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

STEPS="${STEPS:-}"
STOP_ON_FAIL="${STOP_ON_FAIL:-true}"
readonly SKIP_RC=77

STAGE_DIR="${PIPELINE_ROOT}/scripts/10_build_data"
REPORT_TSV="${LOG_DIR}/build_data_report.tsv"
printf 'status\tstep\tdetail\n' > "${REPORT_TSV}"

log "Stage 1 — build data (STOP_ON_FAIL=${STOP_ON_FAIL})"
[[ -n "${STEPS}" ]] && log "Step filter: ${STEPS}"

# Gate: verify external prerequisites (preflight objects, vendored sources, Ensembl
# cache) before any step runs.
bash "${STAGE_DIR}/check_required_inputs.sh" \
  || die "Stage 1 aborted: required inputs missing (see above)."

# The ordered build sequence. Adding/renaming a step = editing this list + the folder.
STEP_DIRS=(
  01_stage_clinical
  02_stage_pheno
  03_build_molecular_signatures
  04_build_set_to_id
  05_stage_metabolomics_cvs
  06_generate_qc_norm
  07_build_human_feature_to_gene
  08_build_qc_objects
  09_build_da
  10_build_da_assemble
  11_build_sum_stats
  12_build_camera_results
  13_build_fcm
  14_build_utoronto_tfs
  15_build_ptmsea_input
  16_build_scion
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
  logfile="${LOG_DIR}/data_${step}.log"
  ran=$((ran+1))
  set +e
  run_step "10_build_data/${step}" "${logfile}" bash "${bs}"
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

log "Stage 1: ${ran} ran, ${passed} PASS, ${skipped} SKIP, ${failed} FAIL — see ${REPORT_TSV}"
(( failed == 0 )) || die "Stage 1 finished with ${failed} failed step(s)."
ok "Stage 1 complete."
