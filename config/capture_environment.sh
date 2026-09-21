#!/usr/bin/env bash
# Environment capture: record the machine, the date, the toolchain, the git
# state of every in-scope repo, and every R package version a run could touch.
#
# Not part of the DAG — run it on demand (`make env`) when you want to refresh
# the provenance record, e.g. before committing regenerated outputs or shipping
# a release. Nothing else depends on it.
#
# Output: docs/ENVIRONMENT.md            (checked in — the documentation)
#         docs/environment/*.tsv|.txt    (machine-readable companions)
#         logs/environment_<stamp>.md    (per-run snapshot, gitignored)
#
# Never fails the pipeline: every probe is guarded, missing tools render as "—".

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../scripts/lib/common.sh"

ENV_DIR="${PIPELINE_ROOT}/docs/environment"
OUT_MD="${PIPELINE_ROOT}/docs/ENVIRONMENT.md"
PKG_TSV="${ENV_DIR}/package_versions.tsv"
SESSION_TXT="${ENV_DIR}/session_info.txt"
mkdir -p "${ENV_DIR}"

log "Capturing run environment → ${OUT_MD}"

# Print the first line of `cmd --version`-style output, or an em dash.
ver() {
  local bin="$1"; shift
  have "${bin}" || { printf '—'; return 0; }
  local out
  out="$("$bin" "$@" 2>&1 | head -n1 || true)"
  printf '%s' "${out:-—}"
}

# git state of a checkout: "<branch> @ <sha> (clean|N files dirty)"
git_state() {
  local dir="$1"
  [[ -d "${dir}/.git" ]] || { printf 'not a git checkout'; return 0; }
  local branch sha dirty
  branch="$(git -C "${dir}" rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?')"
  sha="$(git -C "${dir}" rev-parse --short HEAD 2>/dev/null || echo '?')"
  dirty="$(git -C "${dir}" status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
  if [[ "${dirty}" == "0" ]]; then
    printf '`%s` @ `%s` (clean)' "${branch}" "${sha}"
  else
    printf '`%s` @ `%s` (**%s uncommitted file(s)**)' "${branch}" "${sha}" "${dirty}"
  fi
}

# Package Version: field from a repo DESCRIPTION, or an em dash.
desc_version() {
  local d="$1/DESCRIPTION"
  [[ -f "${d}" ]] || { printf '—'; return 0; }
  awk '/^Version:/ {print $2; exit}' "${d}" 2>/dev/null || printf '—'
}

# ---- 1. R inventory (packages + sessionInfo) -------------------------------
R_OK=false
if have "${RSCRIPT}"; then
  # PIPELINE_ROOT is scanned too: this repo's own stage scripts import packages
  # the in-scope repos never mention (ChIPseeker and txdbmaker in
  # 10_build_data/data-raw/lib/qc_helpers.R, for instance). Without it the
  # inventory under-reports what a run actually needs, and anything built from
  # the inventory — the container, for one — comes out missing them.
  if "${RSCRIPT}" "${PIPELINE_ROOT}/config/capture_environment.R" \
       "${PKG_TSV}" "${SESSION_TXT}" \
       "${DATA_PKG_REPO}" "${ANALYSIS_PKG_REPO}" \
       "${PIPELINE_ROOT}" \
       >/dev/null 2>"${LOG_DIR}/capture_environment_R.log"; then
    R_OK=true
  else
    warn "R package inventory failed — see ${LOG_DIR}/capture_environment_R.log"
  fi
else
  warn "Rscript not on PATH — package versions will be omitted"
fi

# Render a TSV as a markdown table, optionally filtering on the `direct` column.
tsv_to_md() {
  local tsv="$1" only_direct="${2:-no}"
  awk -F'\t' -v only="${only_direct}" '
    NR==1 { print "| Package | Version | Source | Built under |";
            print "|---|---|---|---|"; next }
    { if (only == "yes" && $5 != "yes") next
      printf "| `%s` | %s | %s | %s |\n", $1, $2, $3, $4 }
  ' "${tsv}"
}

# ---- 2. Assemble the document ----------------------------------------------
RUN_DATE="$(date +'%Y-%m-%d %H:%M:%S %Z (%z)')"
RUN_EPOCH="$(date +%s)"
STAMP="$(date +'%Y%m%dT%H%M%S')"

{
cat <<EOF
# Run environment

**Generated automatically — do not edit by hand.** Regenerate on demand with \`make env\`; it is not
part of the pipeline DAG and does not re-run on every build. Refresh it when the environment changes
or before committing regenerated outputs: it is the record of *which* machine, *when*, and *with
which package versions* those outputs were produced.

Companion machine-readable files live in \`docs/environment/\`:
\`package_versions.tsv\` (full package inventory) and \`session_info.txt\` (raw R \`sessionInfo()\`).

## Run metadata

| Field | Value |
|---|---|
| Date of run | **${RUN_DATE}** |
| Unix timestamp | ${RUN_EPOCH} |
| User | $(id -un 2>/dev/null || echo '—') |
| Host | $(hostname 2>/dev/null || echo '—') |
| Pipeline root | \`${PIPELINE_ROOT}\` |
| Data versions | current \`v${CURRENT_VERSION}\` → new \`v${NEW_VERSION}\` |
| Locale | ${LANG:-—} |

## Operating system

| Field | Value |
|---|---|
EOF

case "$(uname -s)" in
  Darwin)
    printf '| OS | macOS %s (build %s) |\n' \
      "$(sw_vers -productVersion 2>/dev/null || echo '?')" \
      "$(sw_vers -buildVersion 2>/dev/null || echo '?')"
    printf '| Kernel | %s |\n' "$(uname -sr 2>/dev/null || echo '—')"
    printf '| Architecture | %s |\n' "$(uname -m 2>/dev/null || echo '—')"
    printf '| CPU | %s |\n' "$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo '—')"
    printf '| Cores | %s |\n' "$(sysctl -n hw.ncpu 2>/dev/null || echo '—')"
    printf '| Memory | %s GB |\n' \
      "$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 ))"
    ;;
  Linux)
    printf '| OS | %s |\n' \
      "$( (. /etc/os-release 2>/dev/null && echo "${PRETTY_NAME}") || echo '?')"
    printf '| Kernel | %s |\n' "$(uname -sr 2>/dev/null || echo '—')"
    printf '| Architecture | %s |\n' "$(uname -m 2>/dev/null || echo '—')"
    printf '| CPU | %s |\n' \
      "$(awk -F': ' '/model name/ {print $2; exit}' /proc/cpuinfo 2>/dev/null || echo '—')"
    printf '| Cores | %s |\n' "$(nproc 2>/dev/null || echo '—')"
    printf '| Memory | %s GB |\n' \
      "$(awk '/MemTotal/ {printf "%.0f", $2/1048576}' /proc/meminfo 2>/dev/null || echo '—')"
    ;;
  *)
    printf '| OS | %s |\n' "$(uname -a 2>/dev/null || echo 'unknown')"
    ;;
esac

cat <<EOF

## Toolchain

| Tool | Version | Path |
|---|---|---|
| R | $(ver "${RSCRIPT}" -e 'cat(R.version.string)') | \`$(command -v "${RSCRIPT}" 2>/dev/null || echo '—')\` |
| bash | ${BASH_VERSION:-—} | \`$(command -v bash 2>/dev/null || echo '—')\` |
| make | $(ver make --version) | \`$(command -v make 2>/dev/null || echo '—')\` |
| git | $(ver git --version) | \`$(command -v git 2>/dev/null || echo '—')\` |
| gsutil | $(ver "${GSUTIL}" version) | \`$(command -v "${GSUTIL}" 2>/dev/null || echo '—')\` |
| gcloud | $(ver "${GCLOUD}" version) | \`$(command -v "${GCLOUD}" 2>/dev/null || echo '—')\` |
| curl | $(ver curl --version) | \`$(command -v curl 2>/dev/null || echo '—')\` |

## Repository state

Every repo the pipeline reads or writes, at the commit used for this run.

| Repo | Package version | Git state |
|---|---|---|
| \`motrpac-human-presuspension-repro\` | — | $(git_state "${PIPELINE_ROOT}") |
| \`$(basename "${DATA_PKG_REPO}")\` | $(desc_version "${DATA_PKG_REPO}") | $(git_state "${DATA_PKG_REPO}") |
| \`$(basename "${ANALYSIS_PKG_REPO}")\` | $(desc_version "${ANALYSIS_PKG_REPO}") | $(git_state "${ANALYSIS_PKG_REPO}") |

EOF

if [[ "${R_OK}" == "true" ]]; then
  n_all="$(( $(wc -l < "${PKG_TSV}") - 1 ))"
  n_direct="$(awk -F'\t' 'NR>1 && $5=="yes"' "${PKG_TSV}" | wc -l | tr -d ' ')"

  cat <<EOF
## R package versions

${n_all} packages resolved: ${n_direct} declared directly by the three in-scope repos (or checked by
preflight), the rest pulled in as recursive \`Depends\`/\`Imports\`/\`LinkingTo\`. Versions are what is
installed on this machine right now. Full table, including install library paths:
[\`docs/environment/package_versions.tsv\`](environment/package_versions.tsv).

### Directly declared (${n_direct})

EOF
  tsv_to_md "${PKG_TSV}" yes

  cat <<EOF

<details>
<summary><b>All ${n_all} packages, including transitive dependencies</b></summary>

EOF
  tsv_to_md "${PKG_TSV}" no

  cat <<EOF

</details>

## R session

\`\`\`
EOF
  cat "${SESSION_TXT}"
  printf '```\n'
else
  cat <<'EOF'
## R package versions

_Not captured — `Rscript` was unavailable or the inventory failed on this run.
See `logs/capture_environment_R.log`._
EOF
fi

printf '\n---\nCaptured by `config/capture_environment.sh` on %s.\n' "${RUN_DATE}"

} > "${OUT_MD}"

# Per-run snapshot so a given output can be traced back to its exact environment.
cp "${OUT_MD}" "${LOG_DIR}/environment_${STAMP}.md" 2>/dev/null || true

ok "Environment recorded: ${OUT_MD} (snapshot: logs/environment_${STAMP}.md)"
