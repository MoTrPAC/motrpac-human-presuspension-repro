#!/usr/bin/env bash
# Step 15 — PTMSEA_INPUT: format the phosphoproteomics DA results as a PTM-SEA input .gct.
#
# One row per confidently localized phosphosite, one column per selected contrast, values =
# z-statistics. Reads the step-10 *_PROT_PH_DA objects and the step-08 *_PROT_PH_QC feature
# metadata (whose confident_site column comes from step 06); preprocess_PTMSEA() is vendored
# alongside — see PTMSEA_INPUT.R for what is adapted.
#
# THE .gct IS AN INPUT, NOT A RESULT. PTM-SEA is run from the Broad Institute's PTM-SEA /
# ssGSEA2.0 Docker image against PTMsigDB, and THAT STEP IS OUT OF SCOPE FOR THIS REPO: this
# build stops at writing a valid .gct. Nothing here pulls or runs the image, nothing
# downstream reads enrichment scores, and there is no reproducibility claim over the
# container, its PTMsigDB release, or its parameters. See README.md for the handoff.
#
# Build parameters, not prompts: PTMSEA_TISSUES (default "muscle adipose"),
# PTMSEA_CONTRAST_TYPE (default "exercise_with_controls"), PTMSEA_CONTRAST_CATEGORY
# (default "EE-CON RE-CON").
#
# Requires steps 06, 08 and 10 — 06 because confident_site is written into the prot-ph
# metadata_features freeze there, and a freeze built before that will fail this step's gate.
# Output: scripts/10_build_data/data/PTMSEA_INPUT_<tissue>.gct + PTMSEA_INPUT.rda
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${RSCRIPT}" "${HERE}/PTMSEA_INPUT.R"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/ptmsea_tests.R"
fi
