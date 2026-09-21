#!/usr/bin/env Rscript
# Tests for UTORONTO_TFs (step 14) — the prot-ph regulator pool.
#
# The expected feature set is re-derived here straight from the vendored extract and
# HUMAN_FEATURE_TO_GENE rather than imported from the builder, so a mistake in the builder's
# filtering or dedup shows up as a failure instead of being mirrored.
#
# The load-bearing assertion is the column set. `subset_qc()` guards its feature filter with
# `if (!is.null(desired_features))`, so an object missing `feature_id` makes
# `.load_scion_matrixes(subset_TFs = TRUE)` silently return every phosphosite instead of
# erroring — the exact state the released .rda is in. That failure is invisible downstream,
# so it is asserted here.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
suppressWarnings(suppressMessages({ library(dplyr) }))
cat("== utoronto tfs tests ==\n")

tf <- load_built("UTORONTO_TFs")
assert("UTORONTO_TFs built", !is.null(tf), "object missing — run step 14")
if (is.null(tf)) finish()

# ── shape / contract ──────────────────────────────────────────────────────────
assert_cols("UTORONTO_TFs", tf, c("feature_id", "gene_symbol"))
assert("exactly two columns", ncol(tf) == 2L,
       paste("got:", paste(colnames(tf), collapse = ", ")))
assert_unique("feature_id is unique", tf$feature_id)
assert_na("feature_id has no NA", tf$feature_id, max_frac = 0)
assert_na("gene_symbol has no NA", tf$gene_symbol, max_frac = 0)
assert("non-empty", nrow(tf) > 0)

# ── independently re-derived expectations ─────────────────────────────────────
src <- file.path(.root, "scripts", "10_build_data", "14_build_utoronto_tfs",
                 "sources", "DatabaseExtract_v_1.01.txt")
if (!file.exists(src)) {
  skip("re-derive from source", "vendored extract absent")
} else {
  raw <- read.csv(src, sep = "\t", check.names = TRUE)[-1]
  assert("source extract has the curated TF call", "Is.TF." %in% colnames(raw))
  tf_genes <- unique(as.character(raw$HGNC.symbol[raw$Is.TF. == "Yes"]))
  report("curated TF genes in the extract", length(tf_genes))

  ftg <- load_built("HUMAN_FEATURE_TO_GENE")
  if (is.null(ftg)) {
    skip("re-derive against HUMAN_FEATURE_TO_GENE", "step 07 object missing")
  } else {
    expected <- ftg %>%
      as.data.frame() %>%
      dplyr::filter(assay == "prot-ph",
                    as.character(gene_symbol) %in% tf_genes) %>%
      dplyr::pull(feature_id) %>%
      as.character() %>%
      unique()
    assert("feature set matches an independent re-derivation",
           setequal(as.character(tf$feature_id), expected),
           sprintf("built %d, expected %d, symmetric diff %d",
                   nrow(tf), length(expected),
                   length(union(setdiff(as.character(tf$feature_id), expected),
                                setdiff(expected, as.character(tf$feature_id))))))

    # every id must be a real prot-ph feature, not a transcript/prot-pr id that leaked in
    ph_ids <- ftg %>% as.data.frame() %>% dplyr::filter(assay == "prot-ph") %>%
      dplyr::pull(feature_id) %>% as.character() %>% unique()
    assert_subset("every feature_id is a prot-ph feature", as.character(tf$feature_id), ph_ids)

    # gene_symbol must be the site's own gene, not an artifact of the many-to-many join
    pairs <- ftg %>% as.data.frame() %>% dplyr::filter(assay == "prot-ph") %>%
      dplyr::mutate(key = paste(as.character(feature_id), as.character(gene_symbol))) %>%
      dplyr::pull(key)
    assert_subset("each (feature_id, gene_symbol) pair is real",
                  paste(as.character(tf$feature_id), as.character(tf$gene_symbol)), pairs)

    assert("every gene_symbol is a curated TF",
           all(as.character(tf$gene_symbol) %in% tf_genes))
  }
}

# ── the regression this object exists to prevent ──────────────────────────────
# The released .rda is the raw 2,765 x 28 download. If step 14 ever produces that shape
# again, subset_TFs goes back to being a no-op.
assert("not the raw download shape", !("Is.TF." %in% colnames(tf)),
       "object still carries the raw extract columns — subset_TFs would be a no-op")

report("prot-ph TF sites", nrow(tf))
report("distinct TF genes represented", dplyr::n_distinct(tf$gene_symbol))

# ── diff vs the shipped package object ────────────────────────────────────────
# Expected to differ: the shipped object is the un-mapped raw extract. Informational.
diff_vs_package("UTORONTO_TFs", tf)

finish()
