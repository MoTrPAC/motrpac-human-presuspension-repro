#!/usr/bin/env bash
# Step 10 — DA_ASSEMBLE: the package-facing differential-analysis objects.
#
# Reads the step-09 DA freeze (staging/freeze/*/da/*.txt) back, attaches the contrast
# metadata from sources/contrast_converter.txt, and reshapes it into CONTRAST_CONVERTER plus
# one {TISSUE}_{OME}_DA object per fitted table (metabolomics platforms stacked per tissue).
# Nothing is refit — see differential_analysis_results.R for what is assembled and for the
# CV dedup that decides which metab platforms survive it.
#
# Requires step 09 to have run. Output: scripts/10_build_data/data/*.rda.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${RSCRIPT}" "${HERE}/differential_analysis_results.R"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/da_assemble_tests.R"
fi
