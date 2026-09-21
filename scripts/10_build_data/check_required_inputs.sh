#!/usr/bin/env bash
# Stage 1 required-inputs gate — run at the START of 10_build_data, before any
# generator. Verifies every pre-existing input the Stage 1 generators consume is
# present, so a run fails fast HERE with one consolidated list rather than a stem
# dying part-way through on a missing file. HARD gate: any missing input aborts.
#
# Scope = inputs that must PRE-EXIST, not objects Stage 1 itself builds:
#   - Stage 0 (preflight) data objects the qc-norm stems load via lib/qc_helpers.R
#   - vendored external raw sources (e.g. the methylcap beta-values from Yongchao Ge,
#     whose raw preprocessing is not part of this repo)
# Add rows to the arrays below as more stems are wired in.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

report_init "${LOG_DIR}/required_inputs_report.tsv"
log "Stage 1 required-inputs check"

# --- Stage 0 (preflight) data objects (label :: absolute path) --------------
declare -a REQUIRED_FILES=(
  "preflight:OME_TISSUE_CODE::${PIPELINE_ROOT}/scripts/00_preflight/data/OME_TISSUE_CODE.rds"
  "preflight:OUTLIERS::${PIPELINE_ROOT}/scripts/00_preflight/data/OUTLIERS.rds"
  "preflight:COVARIATES_FILE::${PIPELINE_ROOT}/scripts/00_preflight/data/COVARIATES_FILE.rds"
  "preflight:QUANTID_BUCKET_FILES::${PIPELINE_ROOT}/scripts/00_preflight/data/QUANTID_BUCKET_FILES.rds"
)
for entry in "${REQUIRED_FILES[@]}"; do
  label="${entry%%::*}"; path="${entry##*::}"
  if [[ -f "${path}" ]]; then
    record_check PASS "${label}" "present"
  else
    record_check FAIL "${label}" "MISSING: ${path#"${PIPELINE_ROOT}"/} (run Stage 0: make preflight data-objects)"
  fi
done

# --- Vendored external raw sources ------------------------------------------
# MethylCap-seq beta-values: external artifact from Yongchao Ge (raw preprocessing
# not in this repo). One file per tissue; consumed only for feature_ids by
# generate_methylcap_qc_norm.R. See sources/methylcap_qc_norm/../README.md.
methyl_dir="${PIPELINE_ROOT}/scripts/10_build_data/06_generate_qc_norm/sources/methylcap_qc_norm"
for tc in t03-edta t06-muscle t11-adipose; do
  if compgen -G "${methyl_dir}/*${tc}*methylcap*beta-values*.txt" >/dev/null; then
    record_check PASS "methylcap:${tc}" "present"
  else
    record_check FAIL "methylcap:${tc}" "MISSING beta-values for ${tc} in sources/methylcap_qc_norm/ (fetch: make sources)"
  fi
done

# MethylCap-seq DA tables: same deal one tier down — fit with a MALAX GLMM in Yongchao Ge's
# external pipeline, so step 09's methylcap stem copies these in rather than modelling them.
# One file per tissue; see sources/README.md in that step.
methyl_da_dir="${PIPELINE_ROOT}/scripts/10_build_data/09_build_da/sources/methylcap_da"
for tc in t03-edta t06-muscle t11-adipose; do
  if compgen -G "${methyl_da_dir}/*${tc}*methylcap*_da_*.txt" >/dev/null; then
    record_check PASS "methylcap_da:${tc}" "present"
  else
    record_check FAIL "methylcap_da:${tc}" "MISSING DA table for ${tc} in 09_build_da/sources/methylcap_da/ (fetch: make sources)"
  fi
done

# ATAC qc-norm + DA: an input only under RERUN_ATAC=FALSE (config/pipeline.env), where the
# step-06 and step-09 stems copy the released tables into the freeze instead of fitting
# them. Under the default TRUE the pipeline produces both itself and neither is required.
if [[ "$(printf '%s' "${RERUN_ATAC:-TRUE}" | tr '[:lower:]' '[:upper:]')" == "FALSE" ]]; then
  atac_qc_dir="${PIPELINE_ROOT}/scripts/10_build_data/06_generate_qc_norm/sources/atac_qc_norm"
  atac_da_dir="${PIPELINE_ROOT}/scripts/10_build_data/09_build_da/sources/atac_da"
  for tc in t05-pbmc t06-muscle; do
    if compgen -G "${atac_qc_dir}/*${tc}*atac*qc-norm*.txt" >/dev/null; then
      record_check PASS "atac:${tc}" "present"
    else
      record_check FAIL "atac:${tc}" "MISSING qc-norm for ${tc} in sources/atac_qc_norm/ (RERUN_ATAC=FALSE; fetch: make sources GROUPS=atac_qc_norm)"
    fi
    if compgen -G "${atac_da_dir}/*${tc}*atac*_da_*.txt" >/dev/null; then
      record_check PASS "atac_da:${tc}" "present"
    else
      record_check FAIL "atac_da:${tc}" "MISSING DA table for ${tc} in 09_build_da/sources/atac_da/ (RERUN_ATAC=FALSE; fetch: make sources GROUPS=atac_da)"
    fi
  done
fi

# METABOLOMICS_CVS: staged from a vendored .rda rather than rebuilt (its upstream generator
# needs a gated bucket read and a step-07 object). See 05_stage_metabolomics_cvs/sources/.
cvs_src="${PIPELINE_ROOT}/scripts/10_build_data/05_stage_metabolomics_cvs/sources/METABOLOMICS_CVS.rda"
if [[ -f "${cvs_src}" ]]; then
  record_check PASS "metabolomics_cvs" "present"
else
  record_check FAIL "metabolomics_cvs" "MISSING: ${cvs_src#"${PIPELINE_ROOT}"/} (see that folder's README)"
fi

# UTORONTO_TFs: the Human TFs database extract (Lambert et al. 2018), a public external
# reference vendored so step 14 builds offline. See 14_build_utoronto_tfs/sources/.
tf_src="${PIPELINE_ROOT}/scripts/10_build_data/14_build_utoronto_tfs/sources/DatabaseExtract_v_1.01.txt"
if [[ -f "${tf_src}" ]]; then
  record_check PASS "utoronto_tfs" "present"
else
  record_check FAIL "utoronto_tfs" "MISSING: ${tf_src#"${PIPELINE_ROOT}"/} (see that folder's README)"
fi

# Local Ensembl v105 annotation cache: TxDb (ChIPseeker) + gene-attribute table,
# built once by sources/ensembl_v105/build_ensembl_v105_cache.R. Consumed by every
# annotation-dependent stem via lib/qc_helpers.R (offline; no live Ensembl).
ens_dir="${PIPELINE_ROOT}/scripts/00_preflight/data-raw/sources/ensembl_v105"
declare -a ENSEMBL_CACHE=(
  "ensembl_v105:txdb::${ens_dir}/txdb_hsapiens_ensembl_v105.sqlite"
  "ensembl_v105:gene-attributes::${ens_dir}/ensembl_v105_gene_attributes.rds"
)
for entry in "${ENSEMBL_CACHE[@]}"; do
  label="${entry%%::*}"; path="${entry##*::}"
  if [[ -f "${path}" ]]; then
    record_check PASS "${label}" "present"
  else
    record_check FAIL "${label}" "MISSING: ${path#"${PIPELINE_ROOT}"/} (build: make sources)"
  fi
done

# Frozen RefMet / KEGG annotation snapshot, built once by
# sources/refmet/build_refmet_cache.R. Consumed by the step 06 metabolomics stem
# via lib/qc_helpers.R (offline; no live Metabolomics Workbench or KEGG query).
# Both artifacts are committed, so absence means a bad checkout, not a stale build.
refmet_dir="${PIPELINE_ROOT}/scripts/00_preflight/data-raw/sources/refmet"
declare -a REFMET_SNAPSHOT=(
  "refmet:name-map::${refmet_dir}/refmet_name_map.txt"
  "refmet:kegg-compounds::${refmet_dir}/kegg_compound_list.txt"
)
for entry in "${REFMET_SNAPSHOT[@]}"; do
  label="${entry%%::*}"; path="${entry##*::}"
  if [[ -f "${path}" ]]; then
    record_check PASS "${label}" "present"
  else
    record_check FAIL "${label}" "MISSING: ${path#"${PIPELINE_ROOT}"/} (rebuild: Rscript ${refmet_dir#"${PIPELINE_ROOT}"/}/build_refmet_cache.R)"
  fi
done

if (( CHECK_FAILS > 0 )); then
  die "${CHECK_FAILS} required input(s) missing — see ${REPORT_TSV}. Fix the above before running Stage 1."
fi
ok "All required Stage 1 inputs present (${CHECK_WARNS} warning(s))."
