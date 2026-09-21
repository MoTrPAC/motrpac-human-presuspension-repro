#!/usr/bin/env bash
# Step 06 — generate qc-norm freeze for every ome.
#
# Runs each stem (a standalone Rscript alongside this file) in dependency order:
# the three stage_* stems, then all generate_*_qc_norm stems, THEN the prot-pr/ph
# imputed stems, which consume the freshly-written prot-pr/ph qc-norm. The stage_*
# stems publish a resource into the freeze rather than computing anything — the
# removed-samples list, the Ensembl v105 TxDb, the RefMet map — so no qc-norm stem
# depends on them and STEMS= may skip them without breaking the rest. Each stem
# reads its raw inputs (quant-id catalog, preflight objects, vendored sources) and
# writes BIC-named freeze files under staging/freeze/; annotation is offline via the
# local Ensembl v105 cache. Runs through EVERY stem even if one fails, then exits
# non-zero if any failed.
#
# Env knobs: STEMS="a b" runs only those (default: all, in order). The stems also
# read RERUN_ATAC (TRUE, the default, refits ATAC here; FALSE stages it from the
# release), VARIANCEPARTITION_PARALLEL_CORES (dream's workers) and
# PARALLEL_CORES (mice's), all from config/pipeline.env.
#
# Output: logs/qcnorm_<stem>.log per stem; logs/qc_norm_report.tsv summary.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"

GEN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Full ordered list: qc-norm stems, then the imputed stems (must follow prot-pr/ph qc-norm).
ALL_STEMS=(
  stage_removed_samples
  stage_ensembl_txdb
  stage_refmet_map
  generate_prot_ol_qc_norm
  generate_prot_pr_qc_norm
  generate_prot_ph_qc_norm
  generate_transcriptomics_qc_norm
  generate_atac_qc_norm
  generate_methylcap_qc_norm
  generate_metab_qc_norm
  generate_clinical_qc_norm
  generate_prot_pr_imputed
  generate_prot_ph_imputed
)
read -r -a STEM_LIST <<< "${STEMS:-${ALL_STEMS[*]}}"

report_init "${LOG_DIR}/qc_norm_report.tsv"
log "Step 06 QC-norm — running ${#STEM_LIST[@]} stem(s)"

# Inputs must all be present before any stem runs (Ensembl cache, preflight objects, etc.).
bash "${PIPELINE_ROOT}/scripts/10_build_data/check_required_inputs.sh" \
  || die "required inputs missing (see above); build them before running QC-norm."

fails=0
for stem in "${STEM_LIST[@]}"; do
  gen="${GEN_DIR}/${stem}.R"
  if [[ ! -f "${gen}" ]]; then
    record_check FAIL "${stem}" "generator not found: ${gen#"${PIPELINE_ROOT}"/}"; fails=$((fails+1)); continue
  fi
  log "run ${stem}"
  logfile="${LOG_DIR}/qcnorm_${stem}.log"
  if "${RSCRIPT}" "${gen}" > "${logfile}" 2>&1; then
    record_check PASS "${stem}" "ok — $(basename "${logfile}")"
  else
    rc=$?
    tail_msg="$(tail -n1 "${logfile}" 2>/dev/null | tr '\t' ' ')"
    record_check FAIL "${stem}" "rc=${rc}: ${tail_msg} (see $(basename "${logfile}"))"; fails=$((fails+1))
  fi
done

log "QC-norm: $(( ${#STEM_LIST[@]} - fails )) PASS, ${fails} FAIL — see ${REPORT_TSV}"
(( fails == 0 )) || die "QC-norm finished with ${fails} failed stem(s)."
ok "QC-norm complete — all ${#STEM_LIST[@]} stems passed."
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${GEN_DIR}/qc_norm_freeze_tests.R"
fi
