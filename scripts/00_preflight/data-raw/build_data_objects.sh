#!/usr/bin/env bash
# Build the leaf (deps=-) data objects from their vendored raw sources.
#
# Each builder writes its object(s) to ../data/ (rds, or gmt.gz for GMTs). This
# driver runs them all and prints a PASS/FAIL/SKIP summary; a single failure does
# not abort the rest. Download-only MSigDB GMTs are vendored as-is (copied).
#
#   Usage:  bash build_data_objects.sh [name ...]     # default: all
#           GSUTIL=gsutil bash build_data_objects.sh
set -uo pipefail

DR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="$(cd "${DR}/.." && pwd)"
DATA="${STAGE}/data"       # .rds data objects ONLY (R-package standard)
OUTPUT="${STAGE}/output"   # everything else: gmt.gz, json
mkdir -p "${DATA}" "${OUTPUT}"
export PREFLIGHT_DR="${DR}" PREFLIGHT_DATA_DIR="${DATA}" PREFLIGHT_OUTPUT_DIR="${OUTPUT}" GSUTIL="${GSUTIL:-gsutil}"
RSCRIPT="${RSCRIPT:-Rscript}"

# builder script (relative to DR) -> label
BUILDERS=(
  "COLORS_ABBREVIATIONS.R"
  "COVARIATES_FILE.R"
  "QUANTID_BUCKET_FILES.R"
  "OME_TISSUE_CODE.R"
  "gmt_processing/GMT_CellMarker.R"
  "gmt_processing/GMT_MitoCarta.R"
  "gmt_processing/GMT_PSP_kinase.R"
  "gmt_processing/GMT_PTMSigDB.R"
  "gmt_processing/GMT_RefMet_subclass.R"
  "OUTLIERS.R"
  "SPLICING_DA.R"
)
# download-only GMTs: vendored prebuilt, copied verbatim (no local rebuild)
PREBUILT_COPY=(
  "c2.all.v2023.2.Hs.symbols.gmt.gz"
  "c5.go.v2023.2.Hs.symbols.gmt.gz"
)

want=("$@")
selected() { [[ ${#want[@]} -eq 0 ]] && return 0; local n; for n in "${want[@]}"; do [[ "$1" == *"$n"* ]] && return 0; done; return 1; }

pass=0; fail=0; skip=0
printf '%-34s %s\n' "OBJECT" "RESULT"
printf '%-34s %s\n' "------" "------"

for rel in "${BUILDERS[@]}"; do
  selected "${rel}" || continue
  log="$(${RSCRIPT} "${DR}/${rel}" 2>&1)"; rc=$?
  last="$(printf '%s\n' "${log}" | tail -n1)"
  if [[ ${rc} -eq 0 ]]; then
    printf '%-34s PASS  %s\n' "${rel##*/}" "${last}"; ((pass++))
  elif printf '%s' "${log}" | grep -qi 'missing source'; then
    printf '%-34s SKIP  %s\n' "${rel##*/}" "${last}"; ((skip++))
  else
    printf '%-34s FAIL  %s\n' "${rel##*/}" "${last}"; ((fail++))
  fi
done

for f in "${PREBUILT_COPY[@]}"; do
  selected "${f}" || { [[ ${#want[@]} -eq 0 ]] || continue; }
  if cp "${DR}/gmt_processing/prebuilt/${f}" "${OUTPUT}/${f}" 2>/dev/null; then
    printf '%-34s PASS  vendored (download-only MSigDB)\n' "${f}"; ((pass++))
  else
    printf '%-34s FAIL  prebuilt missing: %s\n' "${f}" "${f}"; ((fail++))
  fi
done

echo
echo "Built ${pass} PASS, ${skip} SKIP (source absent), ${fail} FAIL"
echo "  rds objects -> ${DATA}"
echo "  gmt.gz/json -> ${OUTPUT}"

# ---- diff built objects vs the version in the package that creates them -----
REPO_ROOT="$(cd "${STAGE}/../.." && pwd)"
GITHUB_ROOT="$(cd "${REPO_ROOT}/.." && pwd)"
export ANALYSIS_PKG_REPO="${ANALYSIS_PKG_REPO:-${GITHUB_ROOT}/MotrpacHumanPreSuspensionAnalysis}"
export DATA_PKG_REPO="${DATA_PKG_REPO:-${GITHUB_ROOT}/MotrpacHumanPreSuspensionData}"
echo
if [[ -d "${ANALYSIS_PKG_REPO}/data" || -d "${DATA_PKG_REPO}/data" ]]; then
  echo "Diff vs source-package data objects:"
  ${RSCRIPT} "${DR}/diff_data_objects.R"
else
  echo "Diff skipped: source package data/ not found (${ANALYSIS_PKG_REPO})"
fi

[[ ${fail} -eq 0 ]]
