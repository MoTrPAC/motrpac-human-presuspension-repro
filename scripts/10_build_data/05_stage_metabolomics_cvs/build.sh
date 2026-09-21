#!/usr/bin/env bash
# Step 05 — stage METABOLOMICS_CVS from the vendored copy in sources/.
# The CVs are calculated by precovid-analyses'
# QC/human-precovid-sed-adu_all_metabolomics_qc-cvs.Rmd; the bucket resources/ .txt, the
# Analysis package's .rda and the vendored copy here are all versions of that Rmd's
# output. The package generator METABOLOMICS_CVS.R only re-downloads the bucket .txt.
# TODO: replace with a local build, once that Rmd is modernized. Staged for now.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE_DATA="${PIPELINE_ROOT}/scripts/10_build_data/data"; mkdir -p "${STAGE_DATA}"

# Vendored once into sources/ rather than re-copied out of the Analysis package on every
# run, so this step does not depend on that checkout being present or on whatever state
# its data/ happens to be in. See sources/README.md to refresh it.
src="${HERE}/sources/METABOLOMICS_CVS.rda"
[[ -f "${src}" ]] || die "missing vendored ${src#"${PIPELINE_ROOT}"/} (see sources/README.md)"
cp "${src}" "${STAGE_DATA}"/
ok "staged METABOLOMICS_CVS -> data/"

# Also emit the table as a freeze resources .txt, which is how the bucket carries it.
"${RSCRIPT}" "${HERE}/METABOLOMICS_CVS_resource.R"

if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/metabolomics_cvs_tests.R"
fi
