#!/usr/bin/env Rscript
# Tests for MOLECULAR_SIGNATURES (step 03). Downstream: run_cameraPR .create_index /
# .prepare_sets subset by database name then unlist sets; default database set is
# setdiff(names, "PTMSIGDB"). Contract: named list of named lists of character vectors.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
cat("== molecular_signatures tests ==\n")

ms <- load_built("MOLECULAR_SIGNATURES")
assert("MOLECULAR_SIGNATURES built", !is.null(ms), "data/MOLECULAR_SIGNATURES.rda missing")
if (!is.null(ms)) {
  assert("is a named list", is.list(ms) && !is.null(names(ms)) && all(nzchar(names(ms))))
  EXPECTED_DB <- c("BIOCARTA", "KEGG_MEDICUS", "PID", "REACTOME", "WP", "GOBP", "GOCC",
                   "GOMF", "MITOCARTA", "PSP", "REFMET", "PTMSIGDB", "CELLMARKER")
  assert("has the 13 expected databases", setequal(names(ms), EXPECTED_DB),
         sprintf("got [%s]", paste(names(ms), collapse = ", ")))
  assert("PTMSIGDB present (cameraPR default drops it)", "PTMSIGDB" %in% names(ms))
  assert("every database is a named list of sets",
         all(vapply(ms, function(db) is.list(db) && length(db) > 0 && !is.null(names(db)), logical(1))))
  assert("every set is a non-empty character vector",
         all(vapply(ms, function(db) all(vapply(db, function(s) is.character(s) && length(s) > 0, logical(1))), logical(1))))
  assert("set names unique within each database",
         all(vapply(ms, function(db) !any(duplicated(names(db))), logical(1))))
}
diff_vs_package("MOLECULAR_SIGNATURES", ms, tolerant = TRUE)
finish()
