#!/usr/bin/env bash
# Step 02 — stage pheno verbatim from the Data package.
# TODO: replace with a local build (pheno.R, depends on cln_objects). Staged for now.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE_DATA="${PIPELINE_ROOT}/scripts/10_build_data/data"; mkdir -p "${STAGE_DATA}"

src="${DATA_PKG_REPO}/data/pheno.rda"
[[ -f "${src}" ]] || die "missing ${src}"
cp "${src}" "${STAGE_DATA}"/
ok "staged pheno -> data/"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/pheno_tests.R"
fi
