#!/usr/bin/env bash
# Container entrypoint: announce which commits this run is reproducing, make sure
# the two mounted source packages are installed, then hand off to the command.
#
# The banner mirrors the "Repository state" table in docs/ENVIRONMENT.md on
# purpose. A container run is only reproducible with respect to the repos bind-
# mounted into it, and those are the one part of the environment the image cannot
# pin — so it says out loud what it got.
set -euo pipefail

GITHUB_ROOT="${GITHUB_ROOT:-/github}"

git_state() {
  local dir="$1"
  # -e, not -d: in a linked worktree .git is a file pointing at the real git dir.
  [[ -e "${dir}/.git" ]] || { printf 'not a git checkout'; return 0; }
  local branch sha dirty
  branch="$(git -C "${dir}" rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?')"
  sha="$(git -C "${dir}" rev-parse --short HEAD 2>/dev/null || echo '?')"
  dirty="$(git -C "${dir}" status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
  if [[ "${dirty}" == "0" ]]; then printf '%s @ %s (clean)' "${branch}" "${sha}"
  else printf '%s @ %s (%s uncommitted)' "${branch}" "${sha}" "${dirty}"; fi
}

# git refuses to read repos owned by another uid; bind mounts from the host are
# exactly that, and every git call above would otherwise fail.
git config --global --add safe.directory '*' 2>/dev/null || true

# Mfuzz imports tcltk, and preflight loads Mfuzz's namespace. Loading tcltk in a
# container with no display fails, so give it a headless one. Cheap, and only
# started if DISPLAY names a local screen that is not already up.
if [[ -n "${DISPLAY:-}" && "${DISPLAY}" == :* ]] && command -v Xvfb >/dev/null 2>&1; then
  screen="${DISPLAY#:}"; screen="${screen%%.*}"
  if [[ ! -e "/tmp/.X11-unix/X${screen}" ]]; then
    Xvfb "${DISPLAY}" -screen 0 1280x1024x24 >/dev/null 2>&1 &
    for _ in 1 2 3 4 5; do
      [[ -e "/tmp/.X11-unix/X${screen}" ]] && break
      sleep 0.3
    done
  fi
fi

echo "─── motrpac-human-presuspension-repro container ──────────────────────────────────────────"
printf '  R           %s\n' "$(R --version 2>/dev/null | head -n1)"
printf '  Bioconductor %s\n' "$(Rscript -e 'cat(as.character(BiocManager::version()))' 2>/dev/null || echo '—')"
printf '  %-34s %s\n' "motrpac-human-presuspension-repro" "$(git_state "${GITHUB_ROOT}/motrpac-human-presuspension-repro")"
# Sibling repos are out of scope by default; report them only when mounted, so
# the banner stays quiet rather than listing NOT MOUNTED lines every run.
for repo in MotrpacHumanPreSuspensionData MotrpacHumanPreSuspensionAnalysis; do
  [[ -d "${GITHUB_ROOT}/${repo}" ]] \
    && printf '  %-34s %s\n' "${repo}" "$(git_state "${GITHUB_ROOT}/${repo}")"
done
echo "───────────────────────────────────────────────────────────────────────"

# Only relevant when the sibling repos have been mounted in; a no-op otherwise.
if [[ "${SKIP_LOCAL_PKG_INSTALL:-0}" != "1" ]]; then
  /usr/local/bin/install_local_pkgs.sh
fi

exec "$@"
