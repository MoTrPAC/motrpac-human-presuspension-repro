#!/usr/bin/env bash
# Step 14 — UTORONTO_TFs: the prot-ph phosphosites whose gene is a curated human
# transcription factor, built from the vendored Human TFs database extract
# (sources/DatabaseExtract_v_1.01.txt) mapped through step 07's HUMAN_FEATURE_TO_GENE.
#
# Two columns, feature_id and gene_symbol. `.load_scion_matrixes()` reads
# UTORONTO_TFs$feature_id to apply subset_TFs = TRUE; the released .rda is the raw 2,765 x 28
# download, which has no feature_id, so subset_TFs is a silent no-op against it rather than
# an error. See UTORONTO_TFs.R.
#
# Requires step 07 (HUMAN_FEATURE_TO_GENE).
# Output: scripts/10_build_data/data/UTORONTO_TFs.rda
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${RSCRIPT}" "${HERE}/UTORONTO_TFs.R"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/utoronto_tests.R"
fi
