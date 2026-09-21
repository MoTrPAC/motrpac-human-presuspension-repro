#!/usr/bin/env bash
# Step 01 — stage the clinical objects (cln_*) verbatim from the Data package.
# TODO: replace with a local build of the gated clinical pipeline
# (clinic_download.R -> clinic.R). Until then these are copied as-is.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE_DATA="${PIPELINE_ROOT}/scripts/10_build_data/data"; mkdir -p "${STAGE_DATA}"

n=$(ls "${DATA_PKG_REPO}"/data/cln_*.rda 2>/dev/null | wc -l | tr -d ' ')
[[ "${n}" -gt 0 ]] || die "no cln_*.rda in ${DATA_PKG_REPO}/data"
cp "${DATA_PKG_REPO}"/data/cln_*.rda "${STAGE_DATA}"/
ok "staged ${n} cln_* objects -> data/"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/clinical_staging_tests.R"
fi
