#!/usr/bin/env Rscript
# Tests for METABOLOMICS_CVS (step 05). Downstream: .eliminate_redundant_metab (load_qc)
# and .prioritize_metab_by_cv_da filter by tissue+assay, join qc_norm on feature_id,
# keep lowest_CV rows, rename to refmet_name. Contract below is what those joins need.
# (The refmet_name ⊆ HUMAN_FEATURE_TO_GENE coupling is checked in the step-07 test,
# where HUMAN_FEATURE_TO_GENE — a downstream object — actually exists.)
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
cat("== metabolomics_cvs tests ==\n")

mc <- load_built("METABOLOMICS_CVS")
assert("METABOLOMICS_CVS built", !is.null(mc), "data/METABOLOMICS_CVS.rda missing")
if (!is.null(mc)) {
  assert("is a data.frame", is.data.frame(mc))
  assert_cols("METABOLOMICS_CVS", mc,
              c("tissue", "assay", "feature_id", "refmet_name", "feature_cv", "lowest_CV"),
              classes = list(feature_cv = "numeric"))
  assert("nrow > 0", nrow(mc) > 0)
  assert("lowest_CV in {yes,no}", all(as.character(mc$lowest_CV) %in% c("yes", "no")),
         sprintf("other values: %s", paste(setdiff(unique(mc$lowest_CV), c("yes","no")), collapse = ", ")))
  assert("tissue in {adipose,blood,muscle}", all(mc$tissue %in% c("adipose", "blood", "muscle")))
  assert("assay is never metab-t-clinical (exempt from CV filtering)",
         !("metab-t-clinical" %in% mc$assay))
  assert("assay values are metab platforms", all(grepl("^metab-[tu]-", mc$assay)))
}
diff_vs_package("METABOLOMICS_CVS", mc, tolerant = TRUE)
finish()
