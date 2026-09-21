#!/usr/bin/env Rscript
# Build the container's R library to match the machine the results were produced
# on, package for package.
#
# The only input is docs/environment/package_versions.tsv, written by `make env`
# (config/capture_environment.R). This script is the executable counterpart of
# that record: whatever versions are installed locally, the image reproduces.
#
# Usage:
#   Rscript install_packages.R --tsv=<package_versions.tsv> --stage=<stage> [--bioc=3.20]
#                              [--report=<build_report.tsv>]
#
# Each stage runs in two passes:
#
#   1. Bulk pass — install the whole set from one repository, letting the resolver
#      sort out dependency order. Fast, and on amd64 largely from binaries.
#   2. Pin pass — for every package whose installed version does not match the
#      recorded one, reinstall that exact version with dependencies = FALSE (the
#      dependency graph is already satisfied by pass 1).
#
# Doing it in that order matters: pinning ~440 packages one at a time in
# dependency order is slow and brittle, while pinning only the ones that actually
# drifted is neither.
#
# Stages, run in this order so each resolves against the library the last left:
#   bioc | cran | github
#
# Exits non-zero if any package with direct == "yes" is missing afterwards.
# Version drift is reported, not fatal — some recorded versions are no longer
# published anywhere, and an unbuildable image helps nobody.

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}

tsv_path <- arg("tsv", "/opt/reference/package_versions.tsv")
stage <- arg("stage", "all")
bioc_version <- arg("bioc", "3.20")
report_path <- arg("report", "/opt/reference/build_report.tsv")

if (!file.exists(tsv_path)) stop("manifest not found: ", tsv_path)
if (!stage %in% c("bioc", "cran", "github", "all")) stop("unknown --stage: ", stage)

manifest <- read.csv(tsv_path, sep = "\t", stringsAsFactors = FALSE)
for (col in c("package", "version", "source", "direct")) {
  if (!col %in% names(manifest)) stop("manifest is missing the '", col, "' column")
}

options(Ncpus = max(1L, parallel::detectCores()))
# download.file's default 60 s timeout aborts the larger GitHub tarballs: PLIER is
# 156 MB and MotrpacRatTraining6moData about 440 MB.
options(timeout = max(3600, getOption("timeout")))
say <- function(...) cat(sprintf("[install:%s] %s\n", stage, paste0(...)))

# capture_environment.R writes versions through as.character(packageVersion(p)),
# which renders every separator as "." — CRAN's "1.4-8" is recorded as "1.4.8".
# Every comparison here goes through package_version(), which treats them alike.
# (remotes::install_version() does the same internally, so passing the recorded
# string straight through resolves the correct tarball.)
# One retry, because a transient network failure mid-pass is indistinguishable from a
# version that does not exist: remotes reports both as "version 'x' is invalid". A whole
# build's worth of pins was lost to a DNS outage that way.
with_retry <- function(what, expr) {
  for (attempt in 1:2) {
    err <- tryCatch({ expr(); NULL }, error = function(e) conditionMessage(e))
    if (is.null(err)) return(invisible(TRUE))
    if (attempt == 1L) {
      say(sprintf("    %s — retrying in 30s", what))
      Sys.sleep(30)
    } else {
      warning(sprintf("%s: %s", what, err), call. = FALSE)
    }
  }
  invisible(FALSE)
}

same_version <- function(a, b) {
  if (is.na(a) || is.na(b)) return(FALSE)
  isTRUE(tryCatch(package_version(a) == package_version(b),
                  error = function(e) identical(a, b)))
}

installed_db <- function() utils::installed.packages()
version_of <- function(p) {
  db <- installed_db()
  if (p %in% rownames(db)) unname(db[p, "Version"]) else NA_character_
}

## ---- Partition the manifest by install source ------------------------------

is_base <- manifest$source == "base"
is_local <- manifest$source == "local source install"
is_github <- grepl("^github:", manifest$source)
is_bioc <- grepl("^Bioconductor", manifest$source)
is_cran <- !(is_base | is_local | is_github | is_bioc)

say(sprintf("manifest: %d packages (%d base, %d local, %d github, %d bioc, %d cran)",
            nrow(manifest), sum(is_base), sum(is_local), sum(is_github),
            sum(is_bioc), sum(is_cran)))

if (any(is_local)) {
  say("deferred to run time (installed from mounted repos, if any): ",
      paste(manifest$package[is_local], collapse = ", "))
}

# rocker builds R --with-recommended-packages, so MASS, Matrix, survival, mgcv and
# friends already sit at exactly the versions R 4.4.1 bundles — which is what the
# manifest records, because the host got them the same way. Leave them alone
# rather than letting a later CRAN upgrade them out from under the reference.
protected <- function() {
  db <- installed_db()
  prio <- db[, "Priority"]
  rownames(db)[!is.na(prio) & prio %in% c("base", "recommended")]
}

## ---- Pin pass ---------------------------------------------------------------

# Returns the rows whose installed version does not match the recorded one.
drifted <- function(rows) {
  keep <- logical(nrow(rows))
  for (i in seq_len(nrow(rows))) {
    got <- version_of(rows$package[[i]])
    keep[[i]] <- !is.na(got) && !same_version(got, rows$version[[i]])
  }
  rows[keep, , drop = FALSE]
}

pin_cran <- function(rows) {
  rows <- rows[!rows$package %in% protected(), , drop = FALSE]
  rows <- drifted(rows)
  if (!nrow(rows)) { say("pin pass: every CRAN version already matches"); return(invisible()) }
  say(sprintf("pin pass: %d CRAN package(s) to move to their recorded versions",
              nrow(rows)))
  for (i in seq_len(nrow(rows))) {
    p <- rows$package[[i]]; want <- rows$version[[i]]
    say(sprintf("  %s %s -> %s", p, version_of(p), want))
    with_retry(sprintf("could not pin %s to %s", p, want),
               function() remotes::install_version(p, version = want, upgrade = "never",
                                                   dependencies = FALSE, quiet = TRUE))
  }
  invisible()
}

# The reference library's Bioconductor packages come from several releases at
# once, so no single release reproduces it. Search a window of releases for the
# one publishing the exact recorded version and install from there.
BIOC_RELEASES <- c("3.19", "3.20", "3.21", "3.22", "3.23")

bioc_repo_urls <- function(rel) {
  sprintf("https://bioconductor.org/packages/%s/%s", rel,
          c("bioc", "data/annotation", "data/experiment"))
}

.bioc_index <- new.env(parent = emptyenv())
bioc_index <- function(url) {
  key <- url
  if (is.null(.bioc_index[[key]])) {
    .bioc_index[[key]] <- tryCatch(utils::available.packages(repos = url),
                                   error = function(e) NULL)
  }
  .bioc_index[[key]]
}

# -> the repo URL publishing `want`, or NA
find_bioc_repo <- function(pkg, want) {
  for (rel in BIOC_RELEASES) {
    for (url in bioc_repo_urls(rel)) {
      ap <- bioc_index(url)
      if (is.null(ap) || !pkg %in% rownames(ap)) next
      if (same_version(ap[pkg, "Version"], want)) return(url)
    }
  }
  NA_character_
}

pin_bioc <- function(rows) {
  rows <- rows[!rows$package %in% protected(), , drop = FALSE]
  rows <- drifted(rows)
  if (!nrow(rows)) { say("pin pass: every Bioconductor version already matches"); return(invisible()) }
  say(sprintf("pin pass: %d Bioconductor package(s) to move to their recorded versions",
              nrow(rows)))
  # A recorded Bioc version compiles against its recorded neighbours, not the
  # snapshot's: ggtree 3.12.0 fails against treeio 1.30 and tidytree 0.4.8. So the CRAN
  # dependencies move first, and the Bioc packages go in dependency order, twice, so a
  # package that failed before its dependency was pinned gets another try.
  apply_holdbacks()
  pin_cran_dependencies(rows$package)
  rows <- rows[match(dependency_order(rows$package), rows$package), , drop = FALSE]
  for (attempt in 1:2) {
    if (attempt == 2L) {
      rows <- drifted(rows)
      if (!nrow(rows)) break
      say(sprintf("pin pass, retry: %d Bioconductor package(s) still off their recorded versions",
                  nrow(rows)))
    }
    for (i in seq_len(nrow(rows))) {
      p <- rows$package[[i]]; want <- rows$version[[i]]
      url <- find_bioc_repo(p, want)
      if (is.na(url)) {
        say(sprintf("  %s %s -> %s NOT PUBLISHED in Bioc %s..%s, leaving as is",
                    p, version_of(p), want, BIOC_RELEASES[[1L]],
                    BIOC_RELEASES[[length(BIOC_RELEASES)]]))
        next
      }
      say(sprintf("  %s %s -> %s (%s)", p, version_of(p), want, basename(dirname(url))))
      with_retry(sprintf("could not pin %s to %s", p, want),
                 function() utils::install.packages(p, repos = url, type = "source",
                                                    dependencies = FALSE))
    }
  }
  invisible()
}

## ---- Build-order holdbacks --------------------------------------------------
#
# Some recorded versions cannot coexist *at install time* even though they coexist
# happily once installed. The reference library was not built in one pass — it
# accumulated — so it contains combinations that cannot be recreated in the
# obvious order.
#
# The known case: ggtree 3.12.0 (and, through it, enrichplot, ChIPseeker and
# clusterProfiler) calls ggplot2:::check_linewidth while byte-compiling. ggplot2
# 4.0 removed that function, so installing ggtree against the recorded ggplot2
# 4.0.2 fails with "object 'check_linewidth' not found". On the host these
# packages work anyway, because they were byte-compiled back when ggplot2 was
# still 3.x and only upgraded around afterwards.
#
# So: install the holdback version first, let the older packages compile against
# it, and let the CRAN pin pass move it to the recorded version afterwards —
# which is exactly the sequence the host went through. The end state matches the
# reference package for package.
#
# 3.5.2, not 3.5.1: the recorded patchwork 1.3.2 imports ggplot2::is_ggplot, first
# exported in 3.5.2, and 3.5.2 still carries check_linewidth. It is the only release
# both sides load against.
HOLDBACKS <- list(
  ggplot2 = list(version = "3.5.2",
                 why = "ggtree 3.12.0 byte-compiles against ggplot2:::check_linewidth, removed in ggplot2 4.0; patchwork 1.3.2 needs is_ggplot, added in 3.5.2"),
  BH = list(version = "1.84.0-0",
            why = "fgsea 1.30.0 and cytolib 2.16.0 compile as C++11, which the Boost in BH 1.87+ no longer supports")
)

apply_holdbacks <- function() {
  db <- installed_db()
  for (p in names(HOLDBACKS)) {
    if (!p %in% manifest$package) next
    recorded <- manifest$version[match(p, manifest$package)]
    hold <- HOLDBACKS[[p]]$version
    if (same_version(recorded, hold)) next          # nothing to hold back
    if (p %in% rownames(db) && same_version(db[p, "Version"], hold)) next
    say(sprintf("holding %s at %s for the build (recorded: %s)", p, hold, recorded))
    say(sprintf("  reason: %s", HOLDBACKS[[p]]$why))
    if (!requireNamespace("remotes", quietly = TRUE)) utils::install.packages("remotes")
    tryCatch(
      remotes::install_version(p, version = hold, upgrade = "never", quiet = TRUE),
      error = function(e)
        warning(sprintf("could not hold %s at %s: %s", p, hold,
                        conditionMessage(e)), call. = FALSE))
  }
}

## ---- Recorded CRAN dependencies ---------------------------------------------

# The hard (Depends/Imports/LinkingTo) dependency closure of `pkgs`, resolved against
# CRAN and the Bioconductor releases that publish the recorded Bioc versions.
.dependency_db <- NULL
dependency_db <- function() {
  if (is.null(.dependency_db)) {
    indices <- c(lapply(unlist(lapply(c("3.19", bioc_version), bioc_repo_urls)), bioc_index),
                 list(tryCatch(utils::available.packages(repos = getOption("repos")),
                               error = function(e) NULL)))
    indices <- Filter(Negate(is.null), indices)
    db <- do.call(rbind, lapply(indices, function(ap)
      ap[, c("Package", "Version", "Depends", "Imports", "LinkingTo"), drop = FALSE]))
    .dependency_db <<- db[!duplicated(db[, "Package"]), , drop = FALSE]
  }
  .dependency_db
}

hard_dependencies <- function(pkgs, recursive = TRUE) {
  tools::package_dependencies(pkgs, db = dependency_db(), recursive = recursive,
                              which = c("Depends", "Imports", "LinkingTo"))
}

dependency_closure <- function(pkgs) unique(unlist(hard_dependencies(pkgs)))

# `pkgs` reordered so each comes after any of the others it depends on.
dependency_order <- function(pkgs) {
  deps <- hard_dependencies(pkgs)
  ordered <- character(0); left <- pkgs
  while (length(left)) {
    ready <- left[vapply(left, function(p) !any(deps[[p]] %in% left), logical(1))]
    if (!length(ready)) ready <- left            # a cycle: keep the given order
    ordered <- c(ordered, ready); left <- setdiff(left, ready)
  }
  ordered
}

# Moves the CRAN packages in that closure to their recorded versions, installing any
# that are missing. Holdbacks and R's own recommended packages are left alone.
pin_cran_dependencies <- function(pkgs) {
  if (!requireNamespace("remotes", quietly = TRUE)) utils::install.packages("remotes")
  deps <- dependency_closure(pkgs)
  rows <- manifest[is_cran & manifest$package %in% deps &
                     !manifest$package %in% c(names(HOLDBACKS), protected()), , drop = FALSE]
  # Retried because a pin can need a dependency pinned later in the same round; stops
  # as soon as a round moves nothing, so an unpublished version is not retried.
  remaining <- Inf
  for (attempt in 1:3) {
    off <- vapply(seq_len(nrow(rows)), function(i)
      !same_version(version_of(rows$package[[i]]), rows$version[[i]]), logical(1))
    todo <- rows[off, , drop = FALSE]
    if (!nrow(todo) || nrow(todo) >= remaining) return(invisible())
    remaining <- nrow(todo)
    say(sprintf("  %d CRAN dependenc%s of the missing packages to their recorded versions",
                nrow(todo), if (nrow(todo) == 1L) "y" else "ies"))
    for (i in seq_len(nrow(todo))) {
      p <- todo$package[[i]]; want <- todo$version[[i]]
      say(sprintf("    %s %s -> %s", p, version_of(p), want))
      tryCatch(
        remotes::install_version(p, version = want, upgrade = "never",
                                 dependencies = NA, quiet = TRUE),
        error = function(e)
          warning(sprintf("could not pin %s to %s: %s", p, want,
                          conditionMessage(e)), call. = FALSE))
    }
  }
  invisible()
}

## ---- Stage: bioc -----------------------------------------------------------

install_bioc <- function() {
  pkgs <- manifest$package[is_bioc]
  missing <- setdiff(pkgs, rownames(installed_db()))
  if (length(missing)) {
    apply_holdbacks()
    say(sprintf("bulk pass: %d Bioconductor packages via release %s",
                length(missing), bioc_version))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager")
    }
    options(BiocManager.check_repositories = FALSE)
    BiocManager::install(version = bioc_version, ask = FALSE, update = FALSE)
    BiocManager::install(missing, ask = FALSE, update = FALSE, version = bioc_version)
    repair_bioc(missing)
  }
  pin_bioc(manifest[is_bioc, , drop = FALSE])
}

# The bulk pass resolves dependencies for the whole set at once, which can pull a
# held-back package back up to its recorded version before a package that needs
# the older one has finished compiling. Anything left missing afterwards gets
# retried one at a time with the holdbacks re-applied and dependencies frozen, so
# nothing shifts underneath it. Repeated until a round installs nothing new,
# which also sorts out ordering (ggtree before enrichplot before clusterProfiler).
#
# Holding ggplot2 back is not enough on its own. The bulk pass takes every CRAN
# dependency from the snapshot's current release, not the recorded one: aplot 0.2.9
# and ggfun 0.2.0 where the reference has 0.2.3 and 0.1.7, and ggraph, shadowtext and
# downloader not at all, since the retry below freezes dependencies. Each round therefore first moves the CRAN packages the missing ones
# depend on to their recorded versions, installing any that are absent, so the
# retry compiles against the same neighbours the reference library did.
repair_bioc <- function(wanted) {
  for (round in 1:5) {
    still <- setdiff(wanted, rownames(installed_db()))
    if (!length(still)) return(invisible())
    say(sprintf("repair pass %d: %d Bioconductor package(s) still missing", round,
                length(still)))
    apply_holdbacks()
    pin_cran_dependencies(still)
    progress <- FALSE
    for (p in still) {
      want <- manifest$version[match(p, manifest$package)]
      url <- find_bioc_repo(p, want)
      if (is.na(url)) url <- BiocManager::repositories()[["BioCsoft"]]
      tryCatch(
        utils::install.packages(p, repos = url, type = "source", dependencies = FALSE),
        error = function(e) warning(sprintf("repair failed for %s: %s", p,
                                            conditionMessage(e)), call. = FALSE))
      if (p %in% rownames(installed_db())) progress <- TRUE
    }
    if (!progress) {
      say("repair pass made no progress; leaving the rest to the reconcile step")
      return(invisible())
    }
  }
  invisible()
}

## ---- Stage: cran -----------------------------------------------------------

install_cran <- function() {
  pkgs <- manifest$package[is_cran]
  missing <- setdiff(pkgs, rownames(installed_db()))
  if (length(missing)) {
    say(sprintf("bulk pass: %d CRAN packages from %s",
                length(missing), paste(getOption("repos"), collapse = " ")))
    utils::install.packages(missing)
  }
  if (!requireNamespace("remotes", quietly = TRUE)) utils::install.packages("remotes")
  # A package since removed from CRAN (inspectdf, which MotrpacBicQC imports) is not
  # in the snapshot's index, so the bulk pass skips it. remotes finds the recorded
  # version in the snapshot's Archive.
  archived <- manifest[is_cran & !manifest$package %in% rownames(installed_db()), , drop = FALSE]
  for (i in seq_len(nrow(archived))) {
    p <- archived$package[[i]]; want <- archived$version[[i]]
    say(sprintf("not in the snapshot index, installing %s %s from the archive", p, want))
    tryCatch(
      remotes::install_version(p, version = want, upgrade = "never",
                               dependencies = NA, quiet = TRUE),
      error = function(e)
        warning(sprintf("could not install %s %s: %s", p, want, conditionMessage(e)),
                call. = FALSE))
  }
  pin_cran(manifest[is_cran, , drop = FALSE])
}

## ---- Stage: github ---------------------------------------------------------

# "github:wgmao/PLIER@fe4e9b2" -> "wgmao/PLIER@fe4e9b2". The commit is the only
# version information these have — they are not on a versioned repository — so it
# is the one pin here that cannot be reconstructed from a version number.
github_ref <- function(source) sub("^github:", "", source)

install_github <- function() {
  rows <- manifest[is_github, , drop = FALSE]
  # Missing, or installed at another version: plotrix arrives from CRAN as a
  # dependency, but the reference has it from a GitHub commit.
  off <- vapply(seq_len(nrow(rows)), function(i)
    !same_version(version_of(rows$package[[i]]), rows$version[[i]]), logical(1))
  rows <- rows[off, , drop = FALSE]
  if (!nrow(rows)) { say("nothing to do"); return(invisible()) }
  if (!requireNamespace("remotes", quietly = TRUE)) utils::install.packages("remotes")

  if (!nzchar(Sys.getenv("GITHUB_PAT"))) {
    say("note: GITHUB_PAT unset — anonymous GitHub API access is rate limited to ",
        "60 requests/hour, which these ", nrow(rows), " installs can exhaust")
  }

  # Two rounds: MotrpacRatTraining6mo sorts before the MotrpacRatTraining6moData it
  # imports, and none of these are in an index dependency_order() could read.
  for (round in 1:2) {
    if (round == 2L) {
      rows <- rows[!rows$package %in% rownames(installed_db()), , drop = FALSE]
      if (!nrow(rows)) break
    }
    for (i in seq_len(nrow(rows))) {
      ref <- github_ref(rows$source[[i]])
      say("install_github(\"", ref, "\")")
      # upgrade="never": the default would silently upgrade the CRAN packages the
      # pin pass just placed at their recorded versions.
      tryCatch(
        remotes::install_github(ref, upgrade = "never", force = TRUE,
                                build_vignettes = FALSE),
        error = function(e) warning("install_github failed for ", ref, ": ",
                                    conditionMessage(e), call. = FALSE))
    }
  }
  invisible()
}

## ---- Run the requested stage(s) --------------------------------------------

if (stage %in% c("bioc", "all")) install_bioc()
if (stage %in% c("cran", "all")) install_cran()
if (stage %in% c("github", "all")) install_github()

## ---- Reconcile, report, and fail loudly on a missing direct dependency -----
#
# install.packages() signals failure with warning(), not stop(): without this
# block a package that failed to compile still yields exit code 0, Docker commits
# the layer, and the break only surfaces at preflight. Reconcile against the
# library rather than trusting the installer.

expected <- manifest[!is_base & !is_local, , drop = FALSE]
if (stage != "all") {
  keep <- switch(stage, bioc = is_bioc, cran = is_cran, github = is_github)
  expected <- manifest[keep, , drop = FALSE]
}

db <- installed_db()
expected$installed_version <- unname(db[match(expected$package, rownames(db)), "Version"])
expected$status <- ifelse(
  is.na(expected$installed_version), "MISSING",
  ifelse(mapply(same_version, expected$installed_version, expected$version),
         "EXACT", "VERSION_DRIFT"))

# Baked into the image: the machine-readable record of how far this build drifted
# from the reference, readable without mounting anything.
if (nzchar(report_path)) {
  dir.create(dirname(report_path), showWarnings = FALSE, recursive = TRUE)
  prior <- if (file.exists(report_path)) {
    read.csv(report_path, sep = "\t", stringsAsFactors = FALSE)
  } else NULL
  out <- rbind(prior[!prior$package %in% expected$package, , drop = FALSE], expected)
  write.table(out[order(tolower(out$package)), ], report_path, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = TRUE)
}

n <- table(factor(expected$status, levels = c("EXACT", "VERSION_DRIFT", "MISSING")))
say(sprintf("%d exact, %d drifted, %d missing (of %d)",
            n[["EXACT"]], n[["VERSION_DRIFT"]], n[["MISSING"]], nrow(expected)))

drift <- expected[expected$status == "VERSION_DRIFT", , drop = FALSE]
if (nrow(drift)) {
  say("version drift (recorded -> installed):")
  cat(paste0("    ", drift$package, " ", drift$version, " -> ", drift$installed_version),
      sep = "\n")
}

missing <- expected[expected$status == "MISSING", , drop = FALSE]
if (!nrow(missing)) quit(status = 0)

missing_direct <- missing[missing$direct == "yes", , drop = FALSE]
missing_indirect <- missing[missing$direct != "yes", , drop = FALSE]

if (nrow(missing_indirect)) {
  say(sprintf("%d transitive package(s) not installed (tolerated — the manifest records the host's closure, not this platform's):",
              nrow(missing_indirect)))
  cat(paste0("    ", missing_indirect$package, " (", missing_indirect$source, ")"),
      sep = "\n")
}

if (nrow(missing_direct)) {
  say(sprintf("FAILED — %d directly-declared package(s) could not be installed:",
              nrow(missing_direct)))
  cat(paste0("    ", missing_direct$package, " ", missing_direct$version,
             " (", missing_direct$source, ")"), sep = "\n")
  quit(status = 1)
}

quit(status = 0)
