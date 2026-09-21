#!/usr/bin/env Rscript
# Inventory the R half of the run environment for docs/ENVIRONMENT.md.
#
# Usage:
#   Rscript capture_environment.R <out_versions.tsv> <out_session_info.txt> [repo ...]
#
# "Packages used" is resolved as: every package declared in the DESCRIPTION of
# each repo passed on the command line, plus every package referenced by
# library()/require()/pkg:: in repos that have no DESCRIPTION,
# plus the load-bearing set preflight checks, expanded to the full recursive
# Depends/Imports/LinkingTo closure over the installed library.

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
}))

if (length(args) < 2L) {
  stop("usage: capture_environment.R <out.tsv> <session_info.txt> [repo ...]")
}
out_tsv <- args[[1L]]
out_si <- args[[2L]]
repos <- if (length(args) > 2L) args[-(1:2)] else character()

`%||%` <- function(x, y) if (is.null(x) || !nzchar(paste(x, collapse = ""))) y else x

# The set preflight treats as load-bearing; kept in sync with R_PKGS in
# 00_preflight.sh so the doc never under-reports what a run actually needs.
CORE <- c("dplyr", "magrittr", "data.table", "tibble", "tidyr", "ggplot2",
          "jsonlite", "here", "devtools", "Biobase", "ComplexHeatmap", "Mfuzz",
          "TMSig", "variancePartition", "MotrpacBicQC")

## ---- Seed the package set -------------------------------------------------

desc_deps <- function(repo) {
  d <- file.path(repo, "DESCRIPTION")
  if (!file.exists(d)) return(character())
  db <- tryCatch(read.dcf(d), error = function(e) NULL)
  if (is.null(db) || !nrow(db)) return(character())
  flds <- intersect(c("Depends", "Imports", "Suggests", "LinkingTo"), colnames(db))
  if (!length(flds)) return(character())
  raw <- unlist(strsplit(paste(stats::na.omit(db[1L, flds]), collapse = ","), ","))
  raw <- trimws(gsub("\\(.*?\\)", "", raw))
  setdiff(raw[nzchar(raw)], c("R", "NA"))
}

# Repos that ship no DESCRIPTION (analysis repos) get scanned instead.
script_deps <- function(repo) {
  files <- list.files(repo, pattern = "[.](R|r|Rmd|rmd|qmd)$",
                      recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("/(renv|packrat|[.]git)/", files)]
  if (!length(files)) return(character())
  txt <- unlist(lapply(files, function(f)
    tryCatch(readLines(f, warn = FALSE), error = function(e) character())))
  if (!length(txt)) return(character())
  grab <- function(pattern) {
    hits <- regmatches(txt, gregexpr(pattern, txt, perl = TRUE))
    hits <- unlist(hits)
    if (!length(hits)) return(character())
    unique(gsub(pattern, "\\1", hits, perl = TRUE))
  }
  found <- unique(c(
    grab("(?:library|require|requireNamespace|loadNamespace)\\(\\s*[\"']?([A-Za-z][A-Za-z0-9._]*)"),
    grab("([A-Za-z][A-Za-z0-9._]*):::?[A-Za-z._]")
  ))
  # `library(x)` / `x::` where x is a loop variable or argument scans as a
  # package name. Drop the usual suspects; the rest is filtered downstream by
  # requiring the name to resolve against the installed library.
  found[nchar(found) > 2L &
        !found %in% c("pkg", "pkgs", "package", "packages", "lib", "lib_name",
                      "libname", "name", "nm", "dep", "deps", "var", "obj")]
}

declared <- CORE   # from a DESCRIPTION or preflight: absence is a real gap
scanned <- character()  # inferred from source: absence may be a false positive
declared_by <- list()
for (repo in repos) {
  pkgs <- desc_deps(repo)
  via <- "DESCRIPTION"
  if (!length(pkgs)) { pkgs <- script_deps(repo); via <- "source scan" }
  if (length(pkgs)) {
    declared_by[[basename(repo)]] <- sprintf("%d (%s)", length(pkgs), via)
    if (via == "DESCRIPTION") declared <- c(declared, pkgs) else scanned <- c(scanned, pkgs)
  }
}
declared <- sort(unique(declared))
scanned <- sort(unique(setdiff(scanned, declared)))
seed <- sort(unique(c(declared, scanned)))

## ---- Expand to the installed recursive closure -----------------------------

db <- utils::installed.packages()
installed <- rownames(db)

closure <- tryCatch(
  unlist(tools::package_dependencies(
    intersect(seed, installed), db = db,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE)),
  error = function(e) character())

all_pkgs <- sort(unique(c(seed, closure)))
missing_declared <- setdiff(declared, installed)
missing_scanned <- setdiff(scanned, installed)
all_pkgs <- intersect(all_pkgs, installed)

## ---- Describe each package -------------------------------------------------

pkg_source <- function(p) {
  d <- tryCatch(utils::packageDescription(p), error = function(e) NULL)
  if (!is.list(d)) return(NA_character_)
  if (identical(d$Priority, "base")) return("base")

  short_sha <- function(x) if (is.null(x)) "" else paste0("@", substr(x, 1L, 7L))
  fields <- paste(c(d$Repository, d$RemoteRepo, d$RemoteUrl, d$biocViews), collapse = " ")

  # Installed straight from GitHub (remotes/devtools) — the sha is the version
  # that matters, since these are not on a versioned repository.
  if (identical(d$RemoteType, "github")) {
    return(paste0("github:", d$RemoteUsername %||% "?", "/",
                  d$RemoteRepo %||% p, short_sha(d$RemoteSha)))
  }

  # Bioconductor, whether resolved through BiocManager or a bioc git remote.
  # Report the release (3.20, 3.22, …) so the doc pins the Bioc/R pairing.
  if (grepl("bioconductor", fields, ignore.case = TRUE) || !is.null(d$biocViews)) {
    rel <- regmatches(fields, regexpr("packages/([0-9]+[.][0-9]+)", fields))
    rel <- sub("packages/", "", rel)
    return(if (length(rel) && nzchar(rel)) paste("Bioconductor", rel) else "Bioconductor")
  }

  if (identical(d$Repository, "CRAN")) return("CRAN")
  if (!is.null(d$Repository) && nzchar(d$Repository)) return(d$Repository)
  if (!is.null(d$RemoteType) && nzchar(d$RemoteType)) {
    return(paste0(d$RemoteType, short_sha(d$RemoteSha)))
  }
  # No repository recorded: installed from a local source tree (the in-scope
  # packages land here). The commit is in the Repository state table above.
  "local source install"
}

info <- data.frame(
  package    = all_pkgs,
  version    = vapply(all_pkgs, function(p) as.character(utils::packageVersion(p)), ""),
  source     = vapply(all_pkgs, pkg_source, ""),
  built_under = unname(db[all_pkgs, "Built"]),
  direct     = ifelse(all_pkgs %in% seed, "yes", "no"),
  library    = unname(db[all_pkgs, "LibPath"]),
  stringsAsFactors = FALSE
)
info <- info[order(tolower(info$package)), ]

write.table(info, out_tsv, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

## ---- Session info ----------------------------------------------------------

bioc <- tryCatch(as.character(BiocManager::version()),
                 error = function(e) "BiocManager not installed")

si <- c(
  paste("R version:      ", R.version.string),
  paste("Platform:       ", R.version$platform),
  paste("Bioconductor:   ", bioc),
  paste("Packages listed:", nrow(info),
        sprintf("(%d directly declared, %d transitive)",
                sum(info$direct == "yes"), sum(info$direct == "no"))),
  if (length(missing_declared))
    paste("DECLARED BUT NOT INSTALLED (a real gap):",
          paste(missing_declared, collapse = ", ")),
  if (length(missing_scanned))
    paste("Referenced in source but not installed (may include false positives",
          "where a variable was scanned as a package name):",
          paste(missing_scanned, collapse = ", ")),
  "",
  "Declared dependencies per repo:",
  if (length(declared_by))
    paste0("  ", names(declared_by), ": ", unlist(declared_by))
  else "  (none resolved)",
  "",
  ".libPaths():",
  paste0("  ", .libPaths()),
  "",
  "sessionInfo():",
  capture.output(sessionInfo())
)
writeLines(si, out_si)
