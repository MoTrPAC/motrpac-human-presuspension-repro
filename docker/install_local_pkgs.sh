#!/usr/bin/env bash
# Install the two in-scope R packages that live in mounted source trees rather
# than in any repository. `docs/environment/package_versions.tsv` records them as
# "local source install", so docker/install_packages.R deliberately skips them at
# build time — they do not exist until docker-compose bind-mounts the sibling
# repos.
#
# They go into ${LOCAL_LIB} (a named volume, not an image layer) so the install
# survives between `docker compose run` invocations: MotrpacHumanPreSuspensionData
# carries ~600 MB of data/, and reinstalling it on every container start would
# dominate the runtime of anything else.
#
# Reinstalls only when the source tree actually changed — the stamp is the
# DESCRIPTION Version plus the git SHA and dirty state. FORCE=1 overrides.
set -euo pipefail

GITHUB_ROOT="${GITHUB_ROOT:-/github}"
LOCAL_LIB="${R_LIBS_USER:-/local-lib}"
FORCE="${FORCE:-0}"

# Analysis before Data: the data package imports the analysis package.
PKGS=(MotrpacHumanPreSuspensionAnalysis MotrpacHumanPreSuspensionData)

mkdir -p "${LOCAL_LIB}"

stamp_of() {
  local repo="$1" version sha dirty
  version="$(awk '/^Version:/ {print $2; exit}' "${repo}/DESCRIPTION" 2>/dev/null || echo '?')"
  if [[ -d "${repo}/.git" ]]; then
    sha="$(git -C "${repo}" rev-parse --short HEAD 2>/dev/null || echo '?')"
    dirty="$(git -C "${repo}" status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
  else
    sha="no-git"; dirty="0"
  fi
  printf '%s %s dirty=%s' "${version}" "${sha}" "${dirty}"
}

for pkg in "${PKGS[@]}"; do
  repo="${GITHUB_ROOT}/${pkg}"
  stamp_file="${LOCAL_LIB}/.installed-${pkg}"

  # Out of scope by default: the container is the environment for motrpac-human-presuspension-repro,
  # and the sibling repos are only mounted when a run needs them.
  if [[ ! -d "${repo}" ]]; then
    continue
  fi

  want="$(stamp_of "${repo}")"
  have="$(cat "${stamp_file}" 2>/dev/null || true)"

  if [[ "${FORCE}" != "1" && "${want}" == "${have}" && -d "${LOCAL_LIB}/${pkg}" ]]; then
    echo "[local-pkgs] ${pkg} up to date (${want})"
    continue
  fi

  echo "[local-pkgs] installing ${pkg} (${want}) -> ${LOCAL_LIB}"
  # Not --no-byte-compile: the analysis package is the one downstream code calls
  # into, and byte-compiling it once here is cheaper than not.
  if R CMD INSTALL --no-docs --library="${LOCAL_LIB}" "${repo}"; then
    printf '%s' "${want}" > "${stamp_file}"
    echo "[local-pkgs] ${pkg} installed"
  else
    rm -f "${stamp_file}"
    echo "[local-pkgs] FAILED to install ${pkg} from ${repo}" >&2
    exit 1
  fi
done
