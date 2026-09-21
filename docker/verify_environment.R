#!/usr/bin/env Rscript
# Compare the R library this script is running against with the reference
# environment recorded in docs/environment/package_versions.tsv.
#
# Usage:
#   Rscript verify_environment.R <reference.tsv> <out.tsv> [<out.md>]
#
# It does not care whether "this library" is inside the container or on a
# workstation, which is the point: `make docker-verify` runs it in the image and
# `make env-diff` runs it natively, and both answer the same question — how does
# what I am about to run with differ from what the released results were produced
# with?
#
# Like config/capture_environment.sh, this never fails the caller. A difference
# is information, not an error.

options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("usage: verify_environment.R <reference.tsv> <out.tsv> [<out.md>]")
}
ref_path <- args[[1L]]
out_tsv <- args[[2L]]
out_md <- if (length(args) > 2L) args[[3L]] else NA_character_

if (!file.exists(ref_path)) stop("reference manifest not found: ", ref_path)
ref <- read.csv(ref_path, sep = "\t", stringsAsFactors = FALSE)

db <- utils::installed.packages()
have <- data.frame(
  package = rownames(db),
  current_version = unname(db[, "Version"]),
  priority = unname(db[, "Priority"]),
  stringsAsFactors = FALSE
)

## ---- Classify --------------------------------------------------------------

cmp <- merge(ref[, c("package", "version", "source", "direct")],
             have[, c("package", "current_version")],
             by = "package", all.x = TRUE)
names(cmp)[names(cmp) == "version"] <- "reference_version"
names(cmp)[names(cmp) == "source"] <- "reference_source"

# Compare as package_version objects, never as strings. capture_environment.R
# writes versions via as.character(packageVersion(p)), which renders every
# separator as "." — so CRAN's "1.4-8" is recorded as "1.4.8". A string compare
# reports 73 spurious differences against the very library the manifest was
# captured from; package_version() treats "." and "-" alike and reports none.
same_version <- function(a, b) {
  mapply(function(x, y) {
    if (is.na(x) || is.na(y)) return(FALSE)
    isTRUE(tryCatch(package_version(x) == package_version(y),
                    error = function(e) x == y))
  }, a, b, USE.NAMES = FALSE)
}

cmp$status <- ifelse(
  is.na(cmp$current_version), "missing",
  ifelse(same_version(cmp$current_version, cmp$reference_version),
         "match", "differs"))

# Packages present here but absent from the reference. Base/recommended are
# excluded: they ship with R and say nothing about how the environment was built.
extra_pkgs <- have[!have$package %in% ref$package &
                     (is.na(have$priority) | have$priority == ""), , drop = FALSE]
extra <- data.frame(
  package = extra_pkgs$package,
  reference_version = NA_character_,
  reference_source = NA_character_,
  direct = "no",
  current_version = extra_pkgs$current_version,
  status = "extra",
  stringsAsFactors = FALSE
)

all_rows <- rbind(cmp[, names(extra)], extra)
all_rows <- all_rows[order(factor(all_rows$status,
                                  levels = c("missing", "differs", "extra", "match")),
                           tolower(all_rows$package)), ]

write.table(all_rows, out_tsv, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

n <- table(factor(all_rows$status, levels = c("match", "differs", "missing", "extra")))
cat(sprintf("[verify] %d match, %d differ, %d missing, %d extra -> %s\n",
            n[["match"]], n[["differs"]], n[["missing"]], n[["extra"]], out_tsv))

missing_direct <- subset(all_rows, status == "missing" & direct == "yes")
if (nrow(missing_direct)) {
  cat(sprintf("[verify] WARNING: %d directly-declared package(s) missing: %s\n",
              nrow(missing_direct), paste(missing_direct$package, collapse = ", ")))
}

if (is.na(out_md)) quit(status = 0)

## ---- Human-readable report -------------------------------------------------

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1L]])) y else x

bioc <- tryCatch(as.character(BiocManager::version()),
                 error = function(e) "not installed")
si <- utils::sessionInfo()
os <- si$running %||% "unknown"
blas <- basename(si$BLAS %||% "unknown")

md_table <- function(df, cols, headers) {
  if (!nrow(df)) return("_None._\n")
  c(paste0("| ", paste(headers, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(headers)), collapse = "|"), "|"),
    apply(df[, cols, drop = FALSE], 1L, function(r)
      paste0("| ", paste(ifelse(is.na(r), "—", r), collapse = " | "), " |")))
}

differs <- subset(all_rows, status == "differs")
differs <- differs[order(differs$direct != "yes", tolower(differs$package)), ]
missing <- subset(all_rows, status == "missing")

lines <- c(
  "# Environment difference report",
  "",
  "**Generated automatically — do not edit by hand.** Regenerate with `make docker-verify`",
  "(inside the container) or `make env-diff` (on this machine). It compares the R library in",
  "use against [`environment/package_versions.tsv`](environment/package_versions.tsv), the",
  "reference environment the released results were produced with.",
  "",
  "## This environment",
  "",
  "| Field | Value |",
  "|---|---|",
  sprintf("| R | %s |", R.version.string),
  sprintf("| Platform | %s |", R.version$platform),
  sprintf("| Running under | %s |", os),
  sprintf("| BLAS | %s |", blas),
  sprintf("| Bioconductor | %s |", bioc),
  sprintf("| Locale | %s |", Sys.getenv("LANG", "—")),
  sprintf("| Time zone | %s |", Sys.timezone()),
  sprintf("| Library paths | %s |", paste(.libPaths(), collapse = ", ")),
  "",
  "## Summary",
  "",
  "| Status | Count | Meaning |",
  "|---|---|---|",
  sprintf("| `match` | %d | same version as the reference |", n[["match"]]),
  sprintf("| `differs` | %d | present, different version |", n[["differs"]]),
  sprintf("| `missing` | %d | in the reference, not installed here |", n[["missing"]]),
  sprintf("| `extra` | %d | installed here, not in the reference |", n[["extra"]]),
  "",
  sprintf("Full row-level detail: [`environment/%s`](environment/%s)",
          basename(out_tsv), basename(out_tsv)),
  "",
  "## Differing versions",
  "",
  "Directly-declared packages are listed first — those are the ones the pipeline names",
  "explicitly, so a version change there is the most likely to change a result.",
  "",
  md_table(differs, c("package", "reference_version", "current_version",
                      "reference_source", "direct"),
           c("Package", "Reference", "Here", "Source", "Direct")),
  "",
  "## Missing",
  "",
  md_table(missing, c("package", "reference_version", "reference_source", "direct"),
           c("Package", "Reference", "Source", "Direct")),
  "",
  "## Differences that are expected by construction",
  "",
  "These will never show as `match`, and chasing them is not useful:",
  "",
  "- **Operating system.** The reference was captured on macOS (arm64, Apple M2 Max). The",
  "  container is Debian bookworm. Anything that reads `Sys.info()` or shells out to a",
  "  platform tool sees a different answer.",
  "- **BLAS/LAPACK.** macOS R links Apple's Accelerate framework; the container uses the",
  "  reference BLAS shipped with R. Floating-point results can differ in the last bits, which",
  "  matters for anything that reports many significant figures or iterates to a tolerance.",
  "- **Bioconductor release.** The reference library is not internally consistent — it spans",
  "  three releases (3.19, 3.20 and 3.21) because it accumulated over time, and its packages",
  "  were built under four different R patch releases (4.4.0 through 4.4.3). The container",
  "  pins a single release, so the packages that came from the other two necessarily shift.",
  "  That is a deliberate improvement, but it is a real difference and it is listed above.",
  "- **Locale and time zone.** The reference ran under `en_US.UTF-8` / `America/Chicago`.",
  "  Sorting of non-ASCII strings and any date formatting follow the locale.",
  "",
  "---",
  sprintf("Generated by `docker/verify_environment.R` on %s.",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)

writeLines(unlist(lines), out_md)
cat(sprintf("[verify] wrote %s\n", out_md))
