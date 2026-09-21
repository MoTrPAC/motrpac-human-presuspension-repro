# Container work — status

**As of 2026-09-17. The image builds to completion and matches the reference: `make docker-verify`
reports 434 match, 0 differ.**

The container is built *from* this repo's own environment capture
(`docs/environment/package_versions.tsv`, written by `make env`), scoped to this repo, pinning the
exact package versions installed on the reference machine.

User-facing documentation lives in the main [`README.md`](../README.md#running-in-docker). This file
is the engineering status — what is proven, what is not, and how to pick it back up.

---

## What exists

| File | Role |
|---|---|
| `Dockerfile` | Ubuntu jammy base (rocker/r-ver:4.4.1), system libraries, Google Cloud SDK, three cacheable R install stages |
| `install_packages.R` | Reads the manifest; per source type runs bulk pass → repair pass → exact-version pin pass |
| `verify_environment.R` | Reference-vs-current package diff; used by both `make docker-verify` and `make env-diff` |
| `entrypoint.sh` | Repo/commit banner, headless X for `tcltk`, optional local package install |
| `install_local_pkgs.sh` | Installs the two mounted source packages into `/local-lib`, stamped so it only reinstalls on change |
| `../docker-compose.yml` | Mounts this repo only; gcloud credentials; writable library volume |
| `../.dockerignore` | Whitelist — keeps ~4.3 GB of vendored data out of the build context |

Plus, outside this folder: `docker-*` and `env-diff` targets in the `Makefile`, the Docker section in
`README.md`, and a bug fix in `config/capture_environment.sh` (see finding 2 below).

---

## Verified

Each of these was checked by running it, not by inspection:

- **Base image facts.** `rocker/r-ver:4.4.1` is Ubuntu 22.04 jammy (not Debian); its default CRAN
  snapshot is `2024-10-30`; on aarch64 it strips the `__linux__/jammy/` path segment.
- **`.dockerignore`** passes exactly the two paths the build needs and nothing else.
- **`make env-diff`** reports 436 match / 0 differ / 0 missing against the library the manifest was
  captured from — i.e. the comparison logic is correct.
- **The whole build, end to end** on `linux/arm64`: all three install stages exit 0 with 72/72
  Bioconductor, 342/342 CRAN and 7/7 GitHub packages at their recorded versions, and
  `make docker-verify` against the built image reports 434 match / 0 differ / 2 missing / 4 extra.
  The 2 missing are the two `local source install` packages, which come from mounted repos at run
  time; the 4 extra are the base image's own (`docopt`, `littler`, `otel`, `ggtangle`).
- **`make docker-preflight`** fails only on `repo:MotrpacHumanPreSuspension{Data,Analysis}` — the
  scope boundary. Tools, R packages, gcloud auth and both buckets pass.
- **CRAN snapshot `2026-07-30`** (the manifest capture date) serves every recorded CRAN version
  except `inspectdf`, which CRAN has since removed and which is installed from the Archive.
- **Bioconductor 3.19–3.23** all resolve; `TMSig` is absent from 3.19 and first appears in 3.20.
- **Bioconductor pin pass** works in both directions — pulled `Biobase` forward from devel (3.23) and
  pushed `ComplexHeatmap` back to 3.19, both onto their recorded versions.
- **`remotes::install_version`** tolerates the manifest's normalized separators — installs `1.4-8`
  from a recorded `1.4.8` (see finding 3).
- **Compose** resolves to this repo only; the parent directory is never mounted.
- **arm64 snapshot branch** of the Dockerfile works in a real build.
- **Entrypoint banner** works against real bind mounts (git's dubious-ownership check is handled).

## Not verified

- Only `linux/arm64` has been built. The `linux/amd64` path is documented but untested end to end.
- No pipeline stage has been run in the container beyond preflight: Stage 1 needs the gated inputs
  and the sibling repos, neither of which the image carries.

---

## Findings worth keeping

These are properties of the reference environment, not of the container. They would bite anyone
handed `docs/ENVIRONMENT.md` and told to match it.

**1. The reference library cannot be rebuilt from scratch in the order its versions imply.**
`ggtree 3.12.0` resolves `ggplot2:::check_linewidth` while byte-compiling; `ggplot2 4.0` removed that
function. Installing the recorded `ggtree` on top of the recorded `ggplot2 4.0.2` fails, taking
`enrichplot`, `ChIPseeker` and `clusterProfiler` with it — two of which are directly declared. On the
reference machine all four load fine, because they were compiled while `ggplot2` was still 3.x and
the upgrade happened around them afterwards. Handled by `HOLDBACKS` in `install_packages.R`: hold
`ggplot2` at 3.5.2 for the Bioconductor stage, then let the CRAN pin pass restore 4.0.2.

3.5.2, not 3.5.1: the recorded `patchwork 1.3.2` imports `ggplot2::is_ggplot`, which 3.5.2 was the
first to export, and 3.5.2 still carries `check_linewidth`. It is the only release both sides load
against.

Holding one package is not enough on its own, which is finding 6.

**2. The environment capture was under-reporting.** `config/capture_environment.sh` passed only the
three sibling repos to the inventory, never `PIPELINE_ROOT`, so packages used solely by this repo's
own stage scripts were invisible — `ChIPseeker` (called at `scripts/10_build_data/data-raw/lib/qc_helpers.R:280`)
and `txdbmaker` among them. Fixed; the manifest went from 428 to 436 packages, 109 to 118 direct.

**3. Versions must be compared as `package_version` objects, never strings.**
`capture_environment.R` writes versions through `as.character(packageVersion(p))`, which renders
every separator as `.` — CRAN's `1.4-8` is recorded as `1.4.8`. A string comparison reported 73
spurious differences against the very library the manifest came from.

**4. The reference library spans multiple Bioconductor releases and four R patch releases.**
Roughly 46 Bioc packages come from 3.19, 22 from 3.20, and a few (`Biobase`, `limma`, `TMSig`) from
devel; packages were built under R 4.4.0 through 4.4.3. No single release pin reproduces it, which is
why the installer searches a release window per package.

**5. Three traps in the base image**, each of which breaks the pipeline silently:
- The `curl` **binary** is purged after R is built (only `libcurl` remains). `00_preflight.sh`
  iterates `for tool in ... curl` and hard-FAILs.
- Tcl/Tk runtime libraries are purged while R still reports `capabilities("tcltk") == TRUE`. `Mfuzz`
  has `Imports: tcltk` and preflight calls `requireNamespace("Mfuzz")`, so Stage 0 would try to load
  Tcl with neither libraries nor a display. Handled with `libtcl8.6`/`libtk8.6` plus an Xvfb display
  started by the entrypoint.
- The default CRAN snapshot predates the manifest by ~21 months.

**6. The bulk pass builds each package against the snapshot's neighbours, not the recorded ones.**
That is a different failure from finding 1, and it bites wherever a recorded version is older than
what the snapshot serves:
- `Deriv 4.3.0` needs an R 4.5 C API (`R_ClosureFormals`) and does not compile on 4.4.1, taking
  `doBy`, `pbkrtest` and `variancePartition` (directly declared) with it. The recorded `Deriv` is
  4.2.0 and builds fine.
- `ggtree 3.12.0` does not byte-compile against `treeio 1.30` / `tidytree 0.4.8`; the recorded
  `treeio` is 1.28.0.
- `fgsea 1.30.0` and `cytolib 2.16.0` compile as C++11, which the Boost in `BH 1.87+` dropped, so
  `BH` is held at 1.84.0-0 alongside `ggplot2` and restored by the CRAN pin pass.

Handled by `pin_cran_dependencies()`: before each repair round and before the Bioconductor pin pass,
the CRAN dependency closure of the packages about to be installed is moved to its recorded versions.
The Bioc pin pass then runs in dependency order, twice.

**7. Three packages cannot be reached the obvious way.**
- `inspectdf` (a hard dependency of `MotrpacBicQC`) is no longer on CRAN. The bulk pass skips what
  the index does not list, so the CRAN stage now installs anything still missing from the Archive by
  its recorded version.
- `PLIER` (156 MB) and `MotrpacRatTraining6moData` (~440 MB) exceed `download.file`'s 60 s default
  timeout; the installer raises it to an hour.
- `plotrix` arrives from CRAN as a transitive dependency, but the reference has it from a GitHub
  commit, so the GitHub stage now also reinstalls packages whose installed version differs from the
  recorded one, and runs twice so `MotrpacRatTraining6mo` can follow its data package.

**8. A network blip is indistinguishable from a version that does not exist.**
`remotes::install_version()` reports both as `version 'x' is invalid for package 'y'`. A DNS outage
mid-build cost 23 pins that way and the build still exited 0, leaving an image that looked fine and
was not. Pins now retry once after 30 s, and `make docker-verify` is the check that matters —
`0 differ` is the only acceptable result.

---

## How to resume

A clean build is around 70 minutes on arm64, so do not iterate on `docker compose build`. Use the
two-step loop:

```bash
# 1. Base image with only the system layers (~40s, cached)
sed -n '1,/^# ---- R packages/p' docker/Dockerfile | sed '$d' > /tmp/base.Dockerfile
docker build -f /tmp/base.Dockerfile -t precovid-base .

# 2. A persistent container to iterate the R install in. Packages already
#    compiled survive between attempts, so a fix costs minutes, not hours.
docker run -d --name precovid-iter \
  -v "$PWD/docs/environment/package_versions.tsv":/opt/reference/package_versions.tsv:ro \
  -v "$PWD/docker/install_packages.R":/usr/local/bin/install_packages.R:ro \
  precovid-base sleep infinity

docker exec precovid-iter bash -c '
  for s in bioc cran github; do
    Rscript /usr/local/bin/install_packages.R --stage=$s --bioc=3.20 \
      --report=/opt/reference/build_report.tsv
  done'
```

Watch for `ERROR: lazy loading failed` / `ERROR: compilation failed` and the `FAILED — N
directly-declared package(s)` line at the end of each stage. Each new install-order conflict gets an
entry in `HOLDBACKS` with its reason. Transitive failures are tolerated by design; only
directly-declared ones fail the build.

Once all three stages run clean in the iteration container:

```bash
make docker-build          # the real image, ~70 min on arm64
make docker-verify         # writes docs/DOCKER_DIFF.md
make docker-preflight      # repo:* FAILs are the scope boundary, not a bug
```

`make docker-verify` is that check: it writes `docs/DOCKER_DIFF.md` and
`docs/environment/docker_diff.tsv` from the image's own library, and `0 differ` is what a good build
looks like. Read it rather than the build log — BuildKit truncates a long step's output, so a stage
summary can be missing from the log entirely.

## Loose ends

- `precovid-base` and the stopped `precovid-iter` container may still exist locally. Remove with
  `docker rm -f precovid-iter && docker rmi precovid-base` if you want a clean slate.
- A rebuild after `make env`, or after editing `install_packages.R`, re-runs all three install
  stages; the tarball cache mount survives, so it is faster than the first build.
- `docs/ENV_DIFF.md` and `docs/environment/env_diff.tsv` are generated by `make env-diff` and are
  currently committed alongside the other generated environment records. Decide whether they belong
  in version control or in `.gitignore`.
- Whether the container should install the two `local source install` packages at all is a scope
  question left open: the mounts are commented out in `docker-compose.yml`, so today it does not.
