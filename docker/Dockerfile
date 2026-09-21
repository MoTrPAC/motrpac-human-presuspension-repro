# motrpac-human-presuspension-repro — the pipeline's software environment, built from its own
# provenance record.
#
# The package set is not hand-maintained here: it is read from
# docs/environment/package_versions.tsv, which `make env`
# (config/capture_environment.R) generates from the machine the released results
# were produced on. Re-run `make env`, rebuild, and the image follows.
#
# What this image is NOT: the data. Individual-level MoTrPAC data, the
# consortium-gated GCS buckets, the vendored data-raw/sources/ trees, and the
# three sibling repos all stay outside and are bind-mounted at run time. See
# docker-compose.yml.
#
# Build:  make docker-build                             (native arch)
#         make docker-build PLATFORM=linux/amd64        (P3M binaries, emulated runtime)
#
# NOTE: rocker/r-ver is Ubuntu, not Debian — 4.x images are built on jammy
# (22.04). Every apt name below is a jammy name.

ARG R_VERSION=4.4.1
FROM rocker/r-ver:${R_VERSION}

# Starting release for the Bioconductor *bulk* pass only — it just needs to pull
# a coherent dependency graph in. The pin pass afterwards searches 3.19 through
# 3.23 for whichever release publishes each package's recorded version, because
# the reference library spans several at once (most of it 3.19 and 3.20, with
# Biobase/TMSig from devel). 3.20 is the best starting point: it matches
# the core S4 stack everything links against, and it is the earliest release
# containing TMSig, which preflight hard-requires.
#
# limma is deliberately held at 3.62.2 (the Bioconductor 3.20 release). It is pinned
# between two hard constraints and 3.62.2 is the only release satisfying both:
#
#   UPPER BOUND — removeBatchEffect() must NOT mean-correct its `covariates` argument.
#   Newer limma subtracts beta %*% t(C - colMeans(C)) instead of beta %*% t(C), which
#   adds the constant beta %*% colMeans(C) to every feature. Five qc-norm stems call it
#   with covariates= (prot-ol, prot-pr, prot-ph, transcriptomics, ATAC), so a newer limma
#   silently re-levels those matrices: structure and all differential contrasts are
#   preserved, absolute values are not.
#
#   Do NOT trust the limma NEWS date here. NEWS attributes this to 3.66.0 (2025-10-28),
#   but 3.64.3 already ships the mean-correction — verified by inspecting the installed
#   function body. The real boundary is between 3.62.2 (clean) and 3.64.3 (corrected).
#
#   LOWER BOUND — limma::eBayes() must accept the `legacy` argument, which
#   variancePartition (1.41.4 here) passes through on every dream() fit. limma 3.60.6
#   predates that argument, so pinning that low breaks every DA stem with
#   "unused argument (legacy = legacy)".
#
# See MotrpacHumanPreSuspensionAnalysis/sandbox/batch_correction_mismatch/limma_version_check/.
# The pin lives in docs/environment/package_versions.tsv; 3.20 is inside the pin pass's
# search window, so recording 3.62.2 there is what enforces it.
ARG BIOC_VERSION=3.20
ENV BIOC_VERSION=${BIOC_VERSION}

# Pin the CRAN snapshot to the manifest capture date (2026-07-30). The base image
# defaults to 2024-10-30, which would silently install e.g. ggplot2 3.5.x where the
# manifest records 4.0.2. Any earlier snapshot cannot serve a recorded version
# released after it: 2026-04-01 lacked cli 3.6.6, rlang, vctrs and six others.
ARG CRAN_SNAPSHOT=2026-07-30

ENV DEBIAN_FRONTEND=noninteractive
# Match the reference run's locale and time zone so string collation and date
# formatting behave the same way. Both are overridable at run time.
ENV LANG=en_US.UTF-8
ENV TZ=America/Chicago

# ---- System libraries -------------------------------------------------------
# Deliberately its own layer with no dependency on the package manifest, so
# regenerating package_versions.tsv does not rebuild any of this.
#
# The bulk of this list is what Posit's sysreqs service returns for exactly the
# CRAN + Bioconductor packages in the manifest. The stragglers are called out
# below, because none of them is guessable from the package list alone.
#
# No LaTeX: every .Rmd in MotrpacPreSuspensionAcute declares output:
# html_document, so pandoc alone covers Stage 3 rendering.
# No libhdf5-dev: Rhdf5lib vendors HDF5, and a system copy risks a header clash.
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates gnupg apt-transport-https \
      curl \
      git make cmake perl python3 \
      pandoc \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      libgit2-dev \
      libicu-dev libsqlite3-dev libuv1-dev \
      zlib1g-dev libbz2-dev liblzma-dev \
      libcairo2-dev libx11-dev \
      libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
      libpng-dev libjpeg-dev libtiff-dev libwebp-dev \
      libglpk-dev libfftw3-dev \
      libgdal-dev gdal-bin libgeos-dev libproj-dev libudunits2-dev libabsl-dev \
      libtcl8.6 libtk8.6 xvfb \
    && rm -rf /var/lib/apt/lists/*
#
# Three of those are not obvious and are load-bearing:
#
#   curl (the binary, not just libcurl4-openssl-dev)
#     rocker purges its build dependencies, and the curl *binary* goes with them.
#     scripts/00_preflight/00_preflight.sh iterates `for tool in ... curl` and
#     hard-FAILs. Without this the container cannot get past Stage 0.
#
#   libtcl8.6 + libtk8.6 + xvfb
#     R is configured --with-tcltk, but tcl-dev/tk-dev are purged and autoremove
#     takes the runtime libraries with them. Mfuzz has `Imports: tcltk`, and
#     preflight calls requireNamespace("Mfuzz"), which loads its imports. So
#     Stage 0 tries to load tcltk in a container that has neither the libraries
#     nor a display. These restore the libraries; DISPLAY below covers the rest.
#
#   libabsl-dev / libwebp-dev / libfftw3-dev / libuv1-dev
#     s2 (via sf), ragg, qqconf and fs respectively. Each is a transitive
#     dependency nothing in the pipeline names directly.

# Headless X display for tcltk. Set as an env var rather than wrapping each
# command in xvfb-run, so `make` targets inside the container need no special
# handling. Xvfb itself is started by the entrypoint.
ENV DISPLAY=:99

# ---- Google Cloud SDK -------------------------------------------------------
# Stage 0 hard-FAILs without gcloud and gsutil on PATH, and every gated input is
# behind a consortium bucket. gsutil is bundled inside google-cloud-cli — there
# is no separate package. Published for both amd64 and arm64.
RUN curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
      | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
      > /etc/apt/sources.list.d/google-cloud-sdk.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends google-cloud-cli \
    && rm -rf /var/lib/apt/lists/*

# ---- CRAN snapshot ----------------------------------------------------------
# rocker bakes its default snapshot into Rprofile.site; append ours so it wins.
# On arm64 the __linux__/<codename>/ segment must be dropped — P3M publishes
# binaries for amd64 only, and asking for the binary path on arm64 yields source
# anyway, just with misleading logs.
RUN . /etc/os-release \
    && if [ "$(uname -m)" = "aarch64" ]; then \
         REPO="https://p3m.dev/cran/${CRAN_SNAPSHOT}"; \
       else \
         REPO="https://p3m.dev/cran/__linux__/${UBUNTU_CODENAME}/${CRAN_SNAPSHOT}"; \
       fi \
    && echo "options(repos = c(CRAN = '${REPO}'), download.file.method = 'libcurl')" \
       >> "${R_HOME}/etc/Rprofile.site" \
    && echo "CRAN snapshot: ${REPO}"

# ---- R packages, from the captured manifest ---------------------------------
COPY docker/install_packages.R /usr/local/bin/install_packages.R
COPY docs/environment/package_versions.tsv /opt/reference/package_versions.tsv

# Three separately-cacheable stages, in dependency order: each resolves against
# the library the previous one left behind. A stage fails the build if any
# directly-declared package in it did not install. The cache mount keeps
# downloaded tarballs across retries, which matters a great deal when a source
# build of ~400 packages fails two hours in.
RUN --mount=type=cache,target=/tmp/Rdownloads,sharing=locked \
    TMPDIR=/tmp/Rdownloads Rscript /usr/local/bin/install_packages.R \
      --tsv=/opt/reference/package_versions.tsv --stage=bioc --bioc=${BIOC_VERSION}
RUN --mount=type=cache,target=/tmp/Rdownloads,sharing=locked \
    TMPDIR=/tmp/Rdownloads Rscript /usr/local/bin/install_packages.R \
      --tsv=/opt/reference/package_versions.tsv --stage=cran --bioc=${BIOC_VERSION}
RUN --mount=type=cache,target=/tmp/Rdownloads,sharing=locked \
    --mount=type=secret,id=github_pat \
    GITHUB_PAT="$(cat /run/secrets/github_pat 2>/dev/null || true)" \
    TMPDIR=/tmp/Rdownloads Rscript /usr/local/bin/install_packages.R \
      --tsv=/opt/reference/package_versions.tsv --stage=github --bioc=${BIOC_VERSION}

# ---- Runtime scripts --------------------------------------------------------
COPY docker/verify_environment.R /usr/local/bin/verify_environment.R
COPY docker/install_local_pkgs.sh /usr/local/bin/install_local_pkgs.sh
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/install_local_pkgs.sh /usr/local/bin/entrypoint.sh

# Bake the image's own difference report in, so `docker run ... cat` can show what
# the image contains without mounting anything.
RUN Rscript /usr/local/bin/verify_environment.R \
      /opt/reference/package_versions.tsv /opt/reference/docker_diff.tsv \
      /opt/reference/DOCKER_DIFF.md || true

# ---- Library layout ---------------------------------------------------------
# Everything above lives in the image's read-only site library. /local-lib is the
# only writable library, backed by a named volume: it holds the two packages
# installed from mounted source trees (install_local_pkgs.sh) and receives
# anything Stage 1's FINALIZE=install writes. Set last so the build above does
# not install into it.
RUN mkdir -p /local-lib
ENV R_LIBS_USER=/local-lib

ENV GITHUB_ROOT=/github
WORKDIR /github/motrpac-human-presuspension-repro

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]
