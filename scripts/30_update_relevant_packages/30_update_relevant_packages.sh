#!/usr/bin/env bash
# Stage 3 — Update the relevant packages.
#
# Carries the regenerated data forward into the two packages that distribute it:
# MotrpacHumanPreSuspensionData (data objects) and
# MotrpacHumanPreSuspensionAnalysis (summary / DA / enrichment objects).
#
# The stage runs in two halves, and only the first is implemented:
#
#   AUDIT (implemented, read-only)  Validate both packages structurally, inventory
#     the code no longer needed, and check that the dependency graph is acyclic.
#     Reads the two package repos, writes only under logs/. A hard circular
#     dependency or a structural ERROR fails the stage at the verdict below, not
#     before the carry — the carry runs either way, and only ever writes under
#     staging/, so a package that cannot install yields a test package that
#     cannot install. AUDIT_ONLY=1 stops here instead of paying for the carry.
#
#   CARRY (implemented, writes only under staging/)  Route every built object to
#     its package, assemble a *test package* per repo from the checkout plus the
#     new data/, re-document, bump versions and NEWS, and run each package's own
#     test suite against the result. Neither package checkout is modified.
#     Promoting a test package into its repo — copying data/, committing, pushing,
#     opening the PR — stays manual, the same posture `make promote` takes toward
#     bucket promotion.
#
# AUDIT_ONLY=1        run the audit and skip the carry entirely.
# STRICT=0            report audit findings but do not fail the stage on them.
# TEST_PKG_ROOT=dir   where the test packages are built (default staging/test-packages).
# BUILD_DATA_DIR=dir  the Stage 1 .rda to carry; defaults to this checkout's, then
#                     to the main checkout's when run from a linked worktree.
# RUN_CHECK=1         also run R CMD check on each test package (slow).
#
# Output: logs/update_relevant_packages_report.tsv
#         logs/package_audit/*.tsv   (per-check detail)
#         logs/package_carry/*.tsv   (routing, manifest, test results)
#         staging/test-packages/<Package>/
# Exit:   non-zero if the audit gate fails under STRICT=1 (the default), or if any
#         carry step fails.

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

AUDIT_DIR_SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/audit"
AUDIT_OUT="${LOG_DIR}/package_audit"
REPORT_TSV="${LOG_DIR}/update_relevant_packages_report.tsv"
STRICT="${STRICT:-1}"
AUDIT_ONLY="${AUDIT_ONLY:-0}"

mkdir -p "${AUDIT_OUT}"
report_init "${REPORT_TSV}"
log "Stage 3 — update relevant packages"

# ---- 1. Are the target packages here? --------------------------------------
# The audit needs both trees on disk. Missing one is fatal to the audit, not just
# worth noting: a half-audited dependency graph cannot answer the cycle question.
PKG_ROOTS=()
for repo in "${DATA_PKG_REPO}" "${ANALYSIS_PKG_REPO}"; do
  if [[ -d "${repo}" ]]; then
    record_check PASS "repo:$(basename "${repo}")" "present at ${repo}"
    PKG_ROOTS+=("${repo}")
  else
    record_check FAIL "repo:$(basename "${repo}")" "not found at ${repo}"
  fi
done

if [[ ${#PKG_ROOTS[@]} -lt 2 ]]; then
  die "both package repos are required for the audit; set GITHUB_ROOT or clone the missing one"
fi

have "${RSCRIPT}" || die "Rscript not on PATH — the audit is R-based"

# Consumer repos: a reference from any of these keeps a function alive. They are
# optional; a missing one shrinks the search scope, so it is recorded, not fatal.
#
# GITHUB_ROOT/MotrpacHumanPreSuspension is deliberately NOT here. Its README
# declares it deprecated and replaced by the Analysis package, and its R/ holds
# older copies of the same functions — so a hit there is a second definition, not
# a call site, and counting it would keep dead code alive forever.
CONSUMER_ARGS=()
for consumer in \
  "${PRECOVID_ROOT}" \
  "${GITHUB_ROOT}/precovid-analyses" \
  "${GITHUB_ROOT}/MotrpacPreSuspensionAcute"; do
  if [[ -d "${consumer}" ]]; then
    CONSUMER_ARGS+=(--consumer "${consumer}")
  else
    record_check WARN "consumer:$(basename "${consumer}")" \
      "not found — call sites there cannot be counted"
  fi
done

# ---- 2. Structure validation -----------------------------------------------
log "auditing package structure -> ${AUDIT_OUT#"${PIPELINE_ROOT}"/}"
if "${RSCRIPT}" "${AUDIT_DIR_SRC}/01_validate_structure.R" \
     "${PKG_ROOTS[@]}" --out "${AUDIT_OUT}" 2>&1 | tee "${LOG_DIR}/audit_structure.log"; then
  record_check PASS "audit:structure" "no ERROR-severity structural findings"
else
  record_check "$([[ "${STRICT}" == "1" ]] && echo FAIL || echo WARN)" \
    "audit:structure" "ERROR-severity findings — see ${AUDIT_OUT#"${PIPELINE_ROOT}"/}/structure_*.tsv"
fi

# ---- 3. Circular dependency check ------------------------------------------
# The gate that matters most. A hard cycle means neither package installs into a
# clean library, so no amount of correct data makes the release usable.
log "checking for circular dependencies"
if "${RSCRIPT}" "${AUDIT_DIR_SRC}/02_check_cycles.R" \
     "${PKG_ROOTS[@]}" --out "${AUDIT_OUT}" 2>&1 | tee "${LOG_DIR}/audit_cycles.log"; then
  record_check PASS "audit:cycles" "no hard circular dependency"
else
  record_check "$([[ "${STRICT}" == "1" ]] && echo FAIL || echo WARN)" \
    "audit:cycles" "circular dependency detected — see ${AUDIT_OUT#"${PIPELINE_ROOT}"/}/cycles_verdict.tsv"
fi

# ---- 4. No-longer-needed code inventory ------------------------------------
# Never a gate. Dead code is a cleanup backlog, not a release blocker, and a
# false DEAD verdict here would otherwise stop a good release.
log "inventorying code no longer needed"
"${RSCRIPT}" "${AUDIT_DIR_SRC}/03_report_unused_code.R" \
  "${PKG_ROOTS[@]}" "${CONSUMER_ARGS[@]}" --out "${AUDIT_OUT}" \
  2>&1 | tee "${LOG_DIR}/audit_unused.log"
record_check PASS "audit:unused-code" \
  "inventory at ${AUDIT_OUT#"${PIPELINE_ROOT}"/}/unused_code_inventory.tsv"

# ---- 5. The carry -----------------------------------------------------------
# Builds a test package per repo under staging/ and carries the objects into it.
# Neither package checkout is written to, here or anywhere else in this stage.
if [[ "${AUDIT_ONLY}" == "1" ]]; then
  log "AUDIT_ONLY=1 — skipping the carry"
else
  CARRY_DIR_SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/carry"
  CARRY_OUT="${LOG_DIR}/package_carry"
  mkdir -p "${CARRY_OUT}"

  # Stage 1 writes its .rda under scripts/10_build_data/data/, which is gitignored
  # — so in a linked worktree that directory is empty while the build sits in the
  # main checkout. Fall back to it rather than failing on an empty source.
  BUILD_DATA_DIR="${BUILD_DATA_DIR:-}"
  if [[ -z "${BUILD_DATA_DIR}" ]]; then
    BUILD_DATA_DIR="${PRECOVID_ROOT}/scripts/10_build_data/data"
    if ! compgen -G "${BUILD_DATA_DIR}/*.rda" >/dev/null; then
      _common="$(git -C "${PRECOVID_ROOT}" rev-parse --path-format=absolute --git-common-dir 2>/dev/null || true)"
      if [[ -n "${_common}" ]]; then
        _main="$(dirname "${_common}")"
        if compgen -G "${_main}/scripts/10_build_data/data/*.rda" >/dev/null; then
          BUILD_DATA_DIR="${_main}/scripts/10_build_data/data"
          log "using the main checkout's Stage 1 build: ${BUILD_DATA_DIR}"
        fi
      fi
    fi
  fi
  # Stage 0's .rds sit beside Stage 1's .rda under the same scripts/ root, so
  # they follow the same main-checkout fallback.
  PREFLIGHT_DATA_DIR="${PREFLIGHT_DATA_DIR:-$(cd "${BUILD_DATA_DIR}/../.." && pwd)/00_preflight/data}"
  TEST_PKG_ROOT="${TEST_PKG_ROOT:-${STAGING_DIR}/test-packages}"
  RUN_CHECK="${RUN_CHECK:-0}"

  carry_step() {
    local n="$1"; shift
    log "carry ${n}"
    if "${RSCRIPT}" "${CARRY_DIR_SRC}/${n}" "$@" 2>&1 | tee "${LOG_DIR}/carry_${n%.R}.log"; then
      record_check PASS "carry:${n%.R}" "ok"
      return 0
    fi
    record_check FAIL "carry:${n%.R}" "see ${LOG_DIR#"${PIPELINE_ROOT}"/}/carry_${n%.R}.log"
    return 1
  }

  carry_step 01_validate_targets.R \
    --data-pkg "${DATA_PKG_REPO}" --analysis-pkg "${ANALYSIS_PKG_REPO}" \
    --build-data "${BUILD_DATA_DIR}" --preflight-data "${PREFLIGHT_DATA_DIR}" \
    --out "${CARRY_OUT}" &&
  carry_step 02_route_objects.R \
    --data-pkg "${DATA_PKG_REPO}" --analysis-pkg "${ANALYSIS_PKG_REPO}" \
    --build-data "${BUILD_DATA_DIR}" --preflight-data "${PREFLIGHT_DATA_DIR}" \
    --inventory "${PIPELINE_ROOT}/docs/data_objects.tsv" --out "${CARRY_OUT}" &&
  carry_step 03_build_test_packages.R \
    --routing "${CARRY_OUT}/routing.tsv" \
    --data-pkg "${DATA_PKG_REPO}" --analysis-pkg "${ANALYSIS_PKG_REPO}" \
    --out-root "${TEST_PKG_ROOT}" --out "${CARRY_OUT}" &&
  carry_step 04_document.R \
    --manifest "${CARRY_OUT}/carry_manifest.tsv" --routing "${CARRY_OUT}/routing.tsv" \
    --data-pkg "${DATA_PKG_REPO}" --analysis-pkg "${ANALYSIS_PKG_REPO}" \
    --out-root "${TEST_PKG_ROOT}" --out "${CARRY_OUT}" &&
  carry_step 05_version_and_news.R \
    --manifest "${CARRY_OUT}/carry_manifest.tsv" --release "${NEW_VERSION}" \
    --data-version "${DATA_PKG_VERSION}" --analysis-version "${ANALYSIS_PKG_VERSION}" \
    --data-pkg "${DATA_PKG_REPO}" --analysis-pkg "${ANALYSIS_PKG_REPO}" \
    --out-root "${TEST_PKG_ROOT}" --out "${CARRY_OUT}" &&
  carry_step 06_check_and_test.R \
    --out-root "${TEST_PKG_ROOT}" --out "${CARRY_OUT}" --check "${RUN_CHECK}"

  log "test packages at ${TEST_PKG_ROOT}"
  warn "Stage 3 carry writes only to staging/ — promoting into the package repos stays manual"
fi

# ---- 6. Verdict -------------------------------------------------------------
if [[ ${CHECK_FAILS} -gt 0 ]]; then
  die "Stage 3: ${CHECK_FAILS} failed check(s), ${CHECK_WARNS} warning(s). See ${REPORT_TSV#"${PIPELINE_ROOT}"/}."
fi
ok "Stage 3 audit complete (${CHECK_WARNS} warning(s)). See ${REPORT_TSV#"${PIPELINE_ROOT}"/}."
