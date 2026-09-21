#!/usr/bin/env bash
# Step 07 — build HUMAN_FEATURE_TO_GENE from the local staging/freeze metadata_features
# (completeness-checked against OME_TISSUE_CODE). Requires step 06 to have run.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${RSCRIPT}" "${HERE}/HUMAN_FEATURE_TO_GENE.R"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/feature_to_gene_tests.R"
fi
