#!/usr/bin/env Rscript
# Tests for the staged clinical objects (step 01). Downstream: load_clinical_data()
# expects each cln_ to be an S3 `motrdat` list with $data (+ $dict). Staged verbatim
# from the Data package, so they must also be identical to the package versions.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
cat("== clinical staging tests ==\n")

cln <- sub("\\.rda$", "", list.files(.T_DATA, pattern = "^cln_.*\\.rda$"))
assert("cln_* objects staged (>= 76)", length(cln) >= 76, sprintf("found %d", length(cln)))

# cln_curated_*/cln_raw_* are S3 `motrdat` lists ($data); cln_chemistry_* are plain
# data.frames (analyte x sample). Accept either.
bad <- character(0)
for (nm in cln) {
  obj <- load_built(nm)
  ok <- !is.null(obj) && ((is.list(obj) && "data" %in% names(obj) && is.data.frame(obj$data)) ||
                          is.data.frame(obj))
  if (!ok) bad <- c(bad, nm)
}
assert("every cln_ is a motrdat list ($data) or a data.frame", length(bad) == 0,
       sprintf("%d bad: %s", length(bad), paste(utils::head(bad, 5), collapse = ", ")))

# staged verbatim -> must equal the Data package version
for (nm in cln) diff_vs_package(nm)
finish()
