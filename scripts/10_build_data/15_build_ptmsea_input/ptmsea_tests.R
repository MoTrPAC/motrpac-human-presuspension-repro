#!/usr/bin/env Rscript
# Tests for PTMSEA_INPUT (step 15) — the PTM-SEA input .gct built from the prot-ph DA
# results. One row per confidently localized phosphosite, one column per selected contrast,
# values = z-statistics.
#
# The expected row set and column set are re-derived here from the inputs (the *_PROT_PH_DA
# objects, the *_PROT_PH_QC feature metadata, and the raw ratio-results confident_site
# column) rather than imported from the builder, so a mistake in the builder's filtering or
# reshaping shows up as a failure instead of being mirrored.
#
# What is NOT tested: anything about PTM-SEA itself. The enrichment is run from the Broad
# Institute's Docker image and is out of scope for this repo, so the contract these tests
# enforce is "a valid, correctly populated .gct exists" and stops there.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
suppressWarnings(suppressMessages({library(dplyr)}))
cat("== ptmsea input tests ==\n")

CONTRAST_TYPE     <- if (nzchar(Sys.getenv("PTMSEA_CONTRAST_TYPE"))) strsplit(trimws(Sys.getenv("PTMSEA_CONTRAST_TYPE")), "[[:space:]]+")[[1]] else "exercise_with_controls"
CONTRAST_CATEGORY <- if (nzchar(Sys.getenv("PTMSEA_CONTRAST_CATEGORY"))) strsplit(trimws(Sys.getenv("PTMSEA_CONTRAST_CATEGORY")), "[[:space:]]+")[[1]] else c("EE-CON", "RE-CON")
TISSUES           <- if (nzchar(Sys.getenv("PTMSEA_TISSUES"))) strsplit(trimws(Sys.getenv("PTMSEA_TISSUES")), "[[:space:]]+")[[1]] else c("muscle", "adipose")

gcts <- load_built("PTMSEA_INPUT")
assert("PTMSEA_INPUT built", !is.null(gcts), "object missing — run step 15")
if (is.null(gcts)) finish()
assert_completeness("every requested tissue was built", names(gcts), TISSUES)

if (!requireNamespace("cmapR", quietly = TRUE)) {
  skip("gct round-trip", "cmapR not installed")
}

# ── independently re-derived expectations ─────────────────────────────────────
# confident_site read straight from the raw ratio-results file, the same source the builder
# uses but resolved here by glob rather than through the quant-id catalog — so a mistake in
# the builder's catalog lookup surfaces as a mismatch rather than being reproduced.
.raw_confident <- function(tissue) {
  tc <- if (tissue == "muscle") "t10-muscle" else "t07-adipose"
  f <- list.files(.T_RAW, pattern = sprintf("%s_prot-ph_ratio-results.*\\.txt$", tc), full.names = TRUE)
  if (length(f) != 1) return(NULL)
  r <- read.csv(f[1], sep = "\t", check.names = FALSE, colClasses = "character",
                stringsAsFactors = FALSE)
  r <- r[!duplicated(r$ptm_id), c("ptm_id", "confident_site"), drop = FALSE]
  stats::setNames(as.logical(toupper(r$confident_site)), r$ptm_id)
}

for (tissue in intersect(TISSUES, names(gcts))) {
  g <- gcts[[tissue]]
  d <- sprintf("%s", tissue)
  assert(sprintf("%s: is a GCT", d), inherits(g, "GCT"))
  if (!inherits(g, "GCT")) next
  m <- g@mat

  assert(sprintf("%s: matrix non-empty", d), nrow(m) > 0 && ncol(m) > 0,
         sprintf("%d x %d", nrow(m), ncol(m)))
  assert(sprintf("%s: matrix is numeric", d), is.numeric(m), sprintf("got %s", class(m[1, 1])))
  assert_unique(sprintf("%s: rids unique", d), g@rid)
  assert(sprintf("%s: rids match matrix rownames", d), identical(g@rid, rownames(m)))
  assert(sprintf("%s: rdesc carries feature_id", d), "feature_id" %in% colnames(g@rdesc))
  # rdesc must be the metadata + the other statistics, never the values themselves.
  assert(sprintf("%s: rdesc has no z.std columns", d),
         !any(grepl("^z\\.std", colnames(g@rdesc))),
         paste(grep("^z\\.std", colnames(g@rdesc), value = TRUE), collapse = ", "))

  # ── columns: one per selected contrast, tissue-prefixed ─────────────────────
  da <- load_built(sprintf("%s_PROT_PH_DA", toupper(tissue)))
  if (is.null(da)) { skip(sprintf("%s: contrast/feature scope", d), "DA object missing"); next }
  da <- as.data.frame(da)
  sel <- da[as.character(da$contrast_type) %in% CONTRAST_TYPE &
              as.character(da$contrast_category) %in% CONTRAST_CATEGORY, , drop = FALSE]
  expected_cid <- paste(tissue, paste0("z.std_", unique(as.character(sel$contrast_short))), sep = ".")
  assert_completeness(sprintf("%s: every selected contrast is a column", d), g@cid, expected_cid)
  assert(sprintf("%s: no out-of-scope columns", d), length(setdiff(g@cid, expected_cid)) == 0,
         paste(utils::head(setdiff(g@cid, expected_cid), 5), collapse = ", "))
  assert(sprintf("%s: cids match matrix colnames", d), identical(g@cid, colnames(m)))

  # ── rows: DA features that are confidently localized, and only those ─────────
  conf <- .raw_confident(tissue)
  if (is.null(conf)) {
    skip(sprintf("%s: confident-site filter", d), "raw ratio-results not in staging/raw-files")
  } else {
    da_feats <- unique(as.character(sel$feature_id))
    expected_rid <- da_feats[conf[da_feats] %in% TRUE]
    assert_completeness(sprintf("%s: every confident DA site is a row", d), g@rid, expected_rid)
    # The regression this guards: filtering with an NA-bearing confident_site keeps the
    # feature as an all-NA row instead of dropping it, which PTM-SEA reads as real data.
    assert(sprintf("%s: no non-confident site kept", d),
           length(setdiff(g@rid, expected_rid)) == 0,
           sprintf("%d extra, e.g. %s", length(setdiff(g@rid, expected_rid)),
                   paste(utils::head(setdiff(g@rid, expected_rid), 5), collapse = ", ")))
    assert(sprintf("%s: not every row is all-NA", d),
           any(rowSums(!is.na(m)) > 0), "matrix is entirely NA")
  }

  # ── the file itself ─────────────────────────────────────────────────────────
  gct_file <- file.path(.T_DATA, sprintf("PTMSEA_INPUT_%s.gct", tissue))
  assert(sprintf("%s: .gct written", d), file.exists(gct_file), gct_file)
  # appenddim = FALSE — a stamped filename would break any handoff script that names the file.
  stamped <- list.files(.T_DATA, pattern = sprintf("^PTMSEA_INPUT_%s_n[0-9]+x[0-9]+\\.gct$", tissue))
  assert(sprintf("%s: filename not dimension-stamped", d), length(stamped) == 0,
         paste(stamped, collapse = ", "))
  if (file.exists(gct_file)) {
    rt <- tryCatch(cmapR::parse_gctx(gct_file), error = function(e) conditionMessage(e))
    if (is.character(rt)) {
      assert(sprintf("%s: .gct re-parses", d), FALSE, rt)
    } else {
      assert(sprintf("%s: .gct re-parses", d), TRUE)
      assert(sprintf("%s: round-trip dims", d), identical(dim(rt@mat), dim(m)),
             sprintf("%s vs %s", paste(dim(rt@mat), collapse = "x"), paste(dim(m), collapse = "x")))
      assert(sprintf("%s: round-trip rids", d), identical(rt@rid, g@rid))
      # cmapR::write_gct rounds to 4 decimal places (its `precision` default, which upstream's
      # own write_gct() call also takes), so the file is exact only to 5e-5. Asserted against
      # that number rather than a loose tolerance: a real corruption moves values much further,
      # and a change to `precision` should fail here and be a deliberate decision.
      dif <- max(abs(unname(rt@mat) - unname(m)), na.rm = TRUE)
      assert(sprintf("%s: round-trip values (max diff %.1e, write precision 4dp)", d, dif),
             is.finite(dif) && dif <= 5e-5, sprintf("max abs diff %g exceeds 4dp rounding", dif))
      # Missing z-scores must survive as missing. A contrast a site was not tested in has to
      # stay NA in the file, not become 0 — PTM-SEA would score a 0 as a real null result.
      assert(sprintf("%s: round-trip NA pattern", d),
             identical(is.na(unname(rt@mat)), is.na(unname(m))))
    }
  }
  report(sprintf("%s: PTM-SEA input", d),
         sprintf("%d sites x %d contrasts -> %s", nrow(m), ncol(m), basename(gct_file)))
}

# Upstream ships no PTMSEA data object — preprocess_PTMSEA() is a package FUNCTION whose
# caller writes the file — so there is nothing to diff against, by design rather than by gap.
report("diff vs package: PTMSEA_INPUT", "no package .rda (upstream ships a function, not an object)")
report("PTM-SEA enrichment", "out of scope — .gct is handed to the Broad PTM-SEA Docker image; this repo stops here")

finish()
