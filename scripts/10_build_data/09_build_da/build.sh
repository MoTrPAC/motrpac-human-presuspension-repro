#!/usr/bin/env bash
# Step 09 — differential analysis (DA) inputs, one stem per ome.
#
# Each stem is a standalone Rscript alongside this file. It reads the step-08 *_QC objects
# via load_qc_local(), fits the acute mixed model with variancePartition::dream through the
# vendored engine alongside it (see README.md), and writes one BIC-named DA table per tissue x ome under
# staging/freeze/<ome-group>/da/.
#
# Runs through EVERY stem even if one fails, then exits non-zero if any failed — the same
# contract as step 06, and for the same reason: dream fits are expensive, so a failure late
# in the list should not discard the earlier results.
#
# Env knobs: STEMS="a b" runs only those stems. Two more reach the stems rather than this
# driver, and both change what lands in the freeze: EBAYES_LEGACY picks limma's
# empirical-Bayes estimator, so it moves every feature's moderated t and p, and RERUN_ATAC
# decides whether ATAC is refitted here (TRUE, the default) or staged from the release
# (FALSE, which needs the vendored tables — see sources/README.md). VARIANCEPARTITION_PARALLEL_CORES
# sets dream's worker count; PARALLEL_CORES drives mice and does nothing in this step.
#
# SCOPE: proteomics, metabolomics, transcriptomics and ATAC are MODELLED here. methylcap is
# fit with a MALAX GLMM in a separate upstream pipeline that is not part of this repo, so its
# stem copies the vendored external tables into the freeze verbatim instead of fitting
# anything — the same treatment the qc-norm step gives its beta-value matrices. Clinical
# chemistry has its own per-feature DA generator upstream (the shared contrast matrix breaks
# on its scattered missingness) and is not adapted yet.
#
# The ATAC stem reads its qc-norm and sample metadata straight from staging/freeze rather
# than through load_qc_local(): step 08 does not package epigen into *_QC.rda, so there is
# no object for it to load. See the header of generate_atac_da.R.
#
# Output: logs/da_<stem>.log per stem; logs/da_report.tsv summary.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"

GEN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ALL_STEMS=(
  generate_prot_ol_da
  generate_prot_pr_da
  generate_prot_ph_da
  generate_clinical_da
  generate_metab_da
  generate_transcriptomics_da
  generate_atac_da
  generate_methylcap_da
)
read -r -a STEM_LIST <<< "${STEMS:-${ALL_STEMS[*]}}"

report_init "${LOG_DIR}/da_report.tsv"
log "Step 09 DA — running ${#STEM_LIST[@]} stem(s)"

# The DA stems consume step-08 output, so gate on the same required inputs as step 06.
bash "${PIPELINE_ROOT}/scripts/10_build_data/check_required_inputs.sh" \
  || die "required inputs missing (see above); build them before running DA."

fails=0
for stem in "${STEM_LIST[@]}"; do
  gen="${GEN_DIR}/${stem}.R"
  if [[ ! -f "${gen}" ]]; then
    record_check FAIL "${stem}" "generator not found: ${gen#"${PIPELINE_ROOT}"/}"; fails=$((fails+1)); continue
  fi
  log "run ${stem}"
  logfile="${LOG_DIR}/da_${stem}.log"
  if "${RSCRIPT}" "${gen}" > "${logfile}" 2>&1; then
    record_check PASS "${stem}" "ok — $(basename "${logfile}")"
  else
    rc=$?
    tail_msg="$(tail -n1 "${logfile}" 2>/dev/null | tr '\t' ' ')"
    record_check FAIL "${stem}" "rc=${rc}: ${tail_msg} (see $(basename "${logfile}"))"; fails=$((fails+1))
  fi
done

log "DA: $(( ${#STEM_LIST[@]} - fails )) PASS, ${fails} FAIL — see ${REPORT_TSV}"
(( fails == 0 )) || die "DA finished with ${fails} failed stem(s)."
ok "DA complete — all ${#STEM_LIST[@]} stems passed."
# Written as an `if` rather than a `&&` chain because under `set -e` a false `&&` tail would
# itself become the script's non-zero exit status when RUN_TESTS=0.
if [[ "${RUN_TESTS:-1}" == "1" && -f "${GEN_DIR}/da_tests.R" ]]; then
  "${RSCRIPT}" "${GEN_DIR}/da_tests.R"
fi
