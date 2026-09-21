#!/usr/bin/env Rscript
# Tests for pheno (step 02). Downstream: load_pheno()/load_qc merge sample metadata on
# vialLabel and rely on the design columns; subset_qc filters on visitcode/Timepoint.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
cat("== pheno tests ==\n")

ph <- load_built("pheno")
assert("pheno built", !is.null(ph), "data/pheno.rda missing")
if (!is.null(ph)) {
  assert("pheno is a motrdat list with $data", is.list(ph) && "data" %in% names(ph) && is.data.frame(ph$data))
  pd <- ph$data
  assert_cols("pheno$data", pd,
              c("vialLabel", "randomGroupCode", "Sex", "Timepoint", "Age_3_groups", "visitcode", "pid"))
  assert("vialLabel unique (one row per sample)", !any(duplicated(pd$vialLabel)))
  assert("visitcode contains ADU_BAS", "ADU_BAS" %in% as.character(pd$visitcode))
  assert("randomGroupCode values valid",
         all(unique(as.character(pd$randomGroupCode)) %in% c("ADUControl", "ADUEndur", "ADUResist", NA)))
  if (is.factor(pd$Timepoint))
    assert_subset("Timepoint levels valid", levels(pd$Timepoint), TIMEPOINT_LEVELS)
}
diff_vs_package("pheno", ph, tolerant = TRUE)
finish()
