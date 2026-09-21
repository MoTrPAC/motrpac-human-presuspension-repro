#!/usr/bin/env Rscript
# Tests for SET_TO_ID (step 04). Downstream: run_cameraPR left-joins results on `set`
# to attach collection/database/set_id/set_short, then adjusts p-values by collection.
# Contract: data.frame; set values must match flattened MOLECULAR_SIGNATURES set names;
# set_id stable + unique.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
cat("== set_to_id tests ==\n")

sti <- load_built("SET_TO_ID")
assert("SET_TO_ID built", !is.null(sti), "data/SET_TO_ID.rda missing")
if (!is.null(sti)) {
  assert("is a data.frame", is.data.frame(sti))
  assert_cols("SET_TO_ID", sti, c("collection", "database", "set_id", "set", "set_short"),
              classes = list(set_id = "character", set = "character", set_short = "character"))
  assert("nrow > 0", nrow(sti) > 0)
  assert_unique("set_id unique per set", sti$set_id)
  assert("set_id zero-padded to constant width", length(unique(nchar(sti$set_id))) == 1)
  # every `set` must resolve against MOLECULAR_SIGNATURES set names (the join key)
  ms <- load_built("MOLECULAR_SIGNATURES")
  if (!is.null(ms)) {
    ms_sets <- unlist(lapply(ms, names), use.names = FALSE)
    assert_subset("set ⊆ flattened MOLECULAR_SIGNATURES set names", sti$set, ms_sets)
  } else skip("set ⊆ MOLECULAR_SIGNATURES set names", "MOLECULAR_SIGNATURES not built")
}
diff_vs_package("SET_TO_ID", sti, tolerant = TRUE)
finish()
