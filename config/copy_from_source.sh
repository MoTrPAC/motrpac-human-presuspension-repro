#!/usr/bin/env bash
# Provide the staged inputs that exceed GitHub's 100 MB per-file limit and are
# therefore gitignored, so a fresh clone can reach Stage 1:
#   - MethylCap beta-value matrices and DA tables, fetched from the production
#     bucket (staged rather than built — they come out of Yongchao Ge's external
#     pipeline, so a clone has no way to produce them)
#   - the Ensembl v105 TxDb, fetched from the bucket's resources/ tier, which step
#     06's stage_ensembl_txdb stem publishes; built from Ensembl only as a fallback
#
# Not part of the DAG — run it on demand (`make sources`) after a fresh clone, or
# when check_required_inputs.sh reports one of these missing.
#
# Reads from PRODUCTION_BUCKET (config/pipeline.env, currently c${CURRENT_VERSION}),
# the same paths each sources/README.md documents. Writes only into sources/
# directories; nothing else is touched.
#
# Already-present files with a matching byte count are skipped, so re-runs are
# cheap and interrupted transfers resume.

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../scripts/lib/common.sh"

usage() {
  cat <<'EOF'
Usage: copy_from_source.sh [--dry-run] [--force] [group ...]

  --dry-run, -n   List what would be fetched and the total transfer size.
  --force,   -f   Re-fetch even when a local file already matches.
  group           Restrict to one or more groups (default: all).

Groups:
  methylcap_qc_norm   3 beta-value matrices  -> 06_generate_qc_norm/sources/
  methylcap_da        3 MALAX GLMM DA tables -> 09_build_da/sources/
  atac_qc_norm        2 log-cpm matrices     -> 06_generate_qc_norm/sources/  (RERUN_ATAC=FALSE only)
  atac_da             2 dream DA tables      -> 09_build_da/sources/          (RERUN_ATAC=FALSE only)
  ensembl_v105        v105 TxDb from resources/, or built from Ensembl if absent there
EOF
}

DRY_RUN=0
FORCE=0
declare -a WANT=()
while (( $# )); do
  case "$1" in
    -n|--dry-run) DRY_RUN=1 ;;
    -f|--force)   FORCE=1 ;;
    -h|--help)    usage; exit 0 ;;
    -*)           usage >&2; die "unknown option: $1" ;;
    *)            WANT+=("$1") ;;
  esac
  shift
done

# group :: remote glob :: local destination
# The methylcap and atac globs each name their ome: the epigenomics qc-norm and da folders
# hold both, and a group must fetch only its own.
declare -a MANIFEST=(
  "methylcap_qc_norm::${PRODUCTION_BUCKET}/epigenomics/qc-norm/*methylcap*beta-values*.txt::${PIPELINE_ROOT}/scripts/10_build_data/06_generate_qc_norm/sources/methylcap_qc_norm"
  "methylcap_da::${PRODUCTION_BUCKET}/epigenomics/da/*methylcap*_da_*.txt::${PIPELINE_ROOT}/scripts/10_build_data/09_build_da/sources/methylcap_da"
  "atac_qc_norm::${PRODUCTION_BUCKET}/epigenomics/qc-norm/*atac*qc-norm*.txt::${PIPELINE_ROOT}/scripts/10_build_data/06_generate_qc_norm/sources/atac_qc_norm"
  "atac_da::${PRODUCTION_BUCKET}/epigenomics/da/*atac*_da_*.txt::${PIPELINE_ROOT}/scripts/10_build_data/09_build_da/sources/atac_da"
)

# ensembl_v105 is built, not fetched, so it has no MANIFEST row.
declare -a ALL_GROUPS=(methylcap_qc_norm methylcap_da atac_qc_norm atac_da ensembl_v105)

# The ATAC groups are an input only under RERUN_ATAC=FALSE; with the default TRUE the
# step-06 and step-09 stems fit those tables themselves and vendoring them would be wrong.
atac_staged() { [[ "$(printf '%s' "${RERUN_ATAC:-TRUE}" | tr '[:lower:]' '[:upper:]')" == "FALSE" ]]; }

file_size() { stat -f%z "$1" 2>/dev/null || stat -c%s "$1" 2>/dev/null || echo 0; }

human() {
  awk -v b="$1" 'BEGIN{
    split("B KB MB GB TB", u, " "); i = 1
    while (b >= 1024 && i < 5) { b /= 1024; i++ }
    printf (i == 1 ? "%d %s" : "%.1f %s"), b, u[i]
  }'
}

wanted() {
  local g="$1" w
  # A bare `make sources` skips the ATAC groups unless RERUN_ATAC=FALSE makes them an
  # input. Naming one explicitly fetches it either way.
  if (( ${#WANT[@]} == 0 )); then
    [[ "${g}" == atac_* ]] && ! atac_staged && return 1
    return 0
  fi
  for w in "${WANT[@]}"; do [[ "${w}" == "${g}" ]] && return 0; done
  return 1
}

# Reject unknown group names rather than silently fetching nothing.
# The count guard is required: bash 3.2 treats "${WANT[@]}" on an empty array as
# an unbound variable under `set -u`.
if (( ${#WANT[@]} > 0 )); then
  for w in "${WANT[@]}"; do
    found=0
    for g in "${ALL_GROUPS[@]}"; do [[ "${g}" == "${w}" ]] && found=1; done
    (( found )) || die "unknown group: ${w} (see --help)"
  done
fi

# The bucket is only needed for the fetch groups.
if wanted methylcap_qc_norm || wanted methylcap_da || wanted atac_qc_norm || wanted atac_da; then
  have "${GSUTIL}" || die "${GSUTIL} not on PATH"
  gcs_can_read "${PRODUCTION_BUCKET}" \
    || die "cannot read ${PRODUCTION_BUCKET} — consortium access required (${GCLOUD} auth login)"
fi

report_init "${LOG_DIR}/copy_from_source_report.tsv"
log "Staged-input fetch from ${PRODUCTION_BUCKET}"

queued_bytes=0
queued_files=0

for entry in "${MANIFEST[@]}"; do
  group="${entry%%::*}"; rest="${entry#*::}"
  glob="${rest%%::*}"; dest="${rest##*::}"
  wanted "${group}" || continue

  # gsutil ls -l emits "<bytes>  <iso-date>  gs://..." plus a trailing TOTAL line.
  # Read loop rather than mapfile: /bin/bash on macOS is 3.2.
  declare -a listing=()
  while IFS= read -r line; do
    [[ -n "${line}" ]] && listing+=("${line}")
  done < <("${GSUTIL}" ls -l "${glob}" 2>/dev/null | grep -v '^TOTAL:' || true)

  if (( ${#listing[@]} == 0 )); then
    record_check FAIL "${group}" "no files match ${glob}"
    continue
  fi

  (( DRY_RUN )) || mkdir -p "${dest}"
  declare -a fetch=()
  declare -a expect=()

  for line in "${listing[@]}"; do
    read -r bytes _ remote <<<"${line}"
    [[ "${bytes}" =~ ^[0-9]+$ && -n "${remote:-}" ]] || continue
    name="$(basename "${remote}")"
    expect+=("${name}::${bytes}")

    if (( ! FORCE )) && [[ -f "${dest}/${name}" ]] \
       && [[ "$(file_size "${dest}/${name}")" == "${bytes}" ]]; then
      continue
    fi
    fetch+=("${remote}")
    queued_bytes=$((queued_bytes + bytes))
    queued_files=$((queued_files + 1))
    log "  queued ${name} ($(human "${bytes}"))"
  done

  if (( ${#fetch[@]} > 0 )) && (( ! DRY_RUN )); then
    log "${group}: fetching ${#fetch[@]} file(s) → ${dest#"${PIPELINE_ROOT}"/}"
    "${GSUTIL}" -m cp "${fetch[@]}" "${dest}/" \
      || warn "${group}: gsutil cp reported an error; per-file status below"
  fi

  # Report final on-disk state, so a partial transfer shows up as FAIL.
  (( ${#expect[@]} > 0 )) || continue
  for e in "${expect[@]}"; do
    name="${e%%::*}"; bytes="${e##*::}"
    if (( DRY_RUN )) && [[ ! -f "${dest}/${name}" ]]; then
      record_check WARN "${group}:${name}" "would fetch ($(human "${bytes}"))"
    elif [[ -f "${dest}/${name}" ]] && [[ "$(file_size "${dest}/${name}")" == "${bytes}" ]]; then
      record_check PASS "${group}:${name}" "$(human "${bytes}")"
    elif [[ -f "${dest}/${name}" ]]; then
      record_check FAIL "${group}:${name}" \
        "size mismatch: local $(human "$(file_size "${dest}/${name}")") vs remote $(human "${bytes}") — re-run with --force"
    else
      record_check FAIL "${group}:${name}" "missing after fetch"
    fi
  done
done

# The Ensembl v105 TxDb is published under the bucket's resources/ tier (step 06's
# stage_ensembl_txdb stem puts it in the freeze), so it downloads like anything else.
# Building it from Ensembl is the fallback for when the bucket does not carry it yet.
if wanted ensembl_v105; then
  ens_dir="${PIPELINE_ROOT}/scripts/00_preflight/data-raw/sources/ensembl_v105"
  txdb="${ens_dir}/txdb_hsapiens_ensembl_v105.sqlite"
  builder="${ens_dir}/build_ensembl_v105_cache.R"

  # The published TxDb carries a _v<version> suffix from v2.0 on (its entry is under
  # resources in config/file_versions.json); releases through v1.4 published it bare.
  # Glob both and take the last match — sorted, that is the highest version — so the
  # bucket decides which copy exists rather than this script pinning a version it may
  # not carry yet. The local cache under sources/ always keeps the bare name: it is an
  # input to this repo, not a release artifact, and every reader names it that way.
  remote_glob="${PRODUCTION_BUCKET}/resources/txdb_hsapiens_ensembl_v105*.sqlite"
  remote_txdb=""

  if [[ -f "${txdb}" ]] && (( ! FORCE )); then
    # Warm checkout: no bucket round-trip needed.
    record_check PASS "ensembl_v105:txdb" "$(human "$(file_size "${txdb}")")"
    ens_done=1
  fi

  # Size and name of the published resource, both empty if the bucket has no copy. The
  # `|| true` matters: gsutil exits non-zero on a missing object, and pipefail
  # would otherwise abort the script here.
  remote_bytes=""
  if (( ! ${ens_done:-0} )) && have "${GSUTIL}"; then
    remote_line="$("${GSUTIL}" ls -l "${remote_glob}" 2>/dev/null \
      | awk '$1 ~ /^[0-9]+$/ {print $1, $NF}' | tail -n 1 || true)"
    remote_bytes="${remote_line%% *}"
    [[ -n "${remote_line}" ]] && remote_txdb="${remote_line##* }"
  fi

  if (( ${ens_done:-0} )); then
    :

  elif (( DRY_RUN )); then
    if [[ -n "${remote_bytes}" ]]; then
      record_check WARN "ensembl_v105:txdb" "would fetch from resources/ ($(human "${remote_bytes}"))"
    else
      record_check WARN "ensembl_v105:txdb" "not on the bucket — would build via ${builder#"${PIPELINE_ROOT}"/}"
    fi

  elif [[ -n "${remote_bytes}" ]]; then
    mkdir -p "${ens_dir}"
    log "ensembl_v105: fetching $(basename "${remote_txdb}") from resources/ ($(human "${remote_bytes}"))"
    # Named destination, not the directory: the published copy may carry a _v<version>
    # suffix that the local cache path does not.
    "${GSUTIL}" -m cp "${remote_txdb}" "${txdb}" \
      || warn "ensembl_v105: gsutil cp reported an error; status below"
    if [[ -f "${txdb}" ]] && [[ "$(file_size "${txdb}")" == "${remote_bytes}" ]]; then
      record_check PASS "ensembl_v105:txdb" "fetched — $(human "${remote_bytes}")"
    else
      record_check FAIL "ensembl_v105:txdb" "fetch incomplete — re-run with --force"
    fi

  elif ! have "${RSCRIPT}"; then
    record_check FAIL "ensembl_v105:txdb" "not on the bucket and ${RSCRIPT} is not on PATH — cannot build"

  else
    warn "ensembl_v105: resources/ has no TxDb yet — building from Ensembl instead (slow, flaky)"
    ens_log="${LOG_DIR}/ensembl_v105_cache.log"
    if FORCE="${FORCE}" "${RSCRIPT}" "${builder}" > "${ens_log}" 2>&1; then
      record_check PASS "ensembl_v105:txdb" "built — $(human "$(file_size "${txdb}")")"
    elif grep -q 'exhausted [0-9]* attempts' "${ens_log}" 2>/dev/null; then
      record_check FAIL "ensembl_v105:txdb" \
        "Ensembl unreachable across every retry — re-run when your connection is better (see $(basename "${ens_log}"))"
    else
      record_check FAIL "ensembl_v105:txdb" \
        "build failed: $(tail -n1 "${ens_log}" 2>/dev/null | tr '\t' ' ') (see $(basename "${ens_log}"))"
    fi
  fi
fi

if (( DRY_RUN )); then
  log "Dry run: ${queued_files} file(s), $(human "${queued_bytes}") would transfer."
  exit 0
fi

if (( CHECK_FAILS > 0 )); then
  die "${CHECK_FAILS} staged input(s) still missing or incomplete — see ${REPORT_TSV}"
fi
ok "Staged inputs present (${CHECK_WARNS} warning(s)) — see ${REPORT_TSV}"
