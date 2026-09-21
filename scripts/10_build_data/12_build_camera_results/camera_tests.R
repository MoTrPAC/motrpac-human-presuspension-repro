#!/usr/bin/env Rscript
# Tests for CAMERA_RESULTS (step 12) — pre-ranked CAMERA enrichment of the step-10 DA
# results. One row per tissue x ome x contrast x set, carrying the two-sample t-statistic,
# its standard-Normal equivalent, and a BH-adjusted p-value.
#
# The expected shape is re-derived from the inputs (the *_DA objects, MOLECULAR_SIGNATURES,
# SET_TO_ID, CONTRAST_CONVERTER) rather than imported from the builder, so a mistake in the
# builder's scope or joins shows up as a failure instead of being mirrored.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
suppressWarnings(suppressMessages({library(dplyr)}))
cat("== camera results tests ==\n")

cam <- load_built("CAMERA_RESULTS")
assert("CAMERA_RESULTS built", !is.null(cam), "object missing — run step 12")
if (is.null(cam)) finish()
cam <- as.data.frame(cam)

COLS <- c("tissue", "assay", "contrast_type", "contrast", "contrast_short",
          "collection", "database", "set_id", "set", "set_short",
          "set_size", "set_size_DB", "size_ratio", "direction",
          "t", "df", "z.std", "p_value", "adj_p_value")
assert_cols("CAMERA_RESULTS", cam, COLS)
assert("CAMERA_RESULTS: column order", identical(colnames(cam), COLS),
       sprintf("got %s", paste(colnames(cam), collapse = " ")))
assert("CAMERA_RESULTS: non-empty", nrow(cam) > 0)

# ── scope: which tissue x ome were tested ─────────────────────────────────────
# run_cameraPR()'s default selected_omes. The clinical DA objects and epigen are deliberately
# not enrichment-tested: clinical chemistry is 9 analytes, and step 10 does not assemble
# epigen at all.
OMES <- c("transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol", "metab")
da_files <- list.files(.T_DATA, pattern = "_DA\\.rda$", full.names = TRUE)
da_files <- da_files[!grepl("SPLICING_DA\\.rda$", da_files)]
da_scope <- dplyr::bind_rows(lapply(da_files, function(f) {
  e <- new.env(); load(f, envir = e); o <- e[[ls(e)[1]]]
  data.frame(tissue = as.character(o$tissue)[1], assay = as.character(o$assay)[1],
             stringsAsFactors = FALSE)
}))
expected <- da_scope[da_scope$assay %in% OMES, ]
got <- unique(cam[, c("tissue", "assay")])
assert_completeness("every DA tissue x ome in scope was tested",
                    paste(got$tissue, got$assay), paste(expected$tissue, expected$assay))
assert("no out-of-scope ome tested", all(as.character(cam$assay) %in% OMES),
       paste(setdiff(unique(as.character(cam$assay)), OMES), collapse = ", "))
# prot-ph is the regression guard for the step-07 flanking_sequence fix: run_cameraPR keys
# phosphosite enrichment on the flanking sequence, so dropping that column silently reduced
# this object to the four non-phospho omes rather than failing.
assert("prot-ph was tested (needs HUMAN_FEATURE_TO_GENE$flanking_sequence)",
       "prot-ph" %in% as.character(cam$assay))

# ── databases ─────────────────────────────────────────────────────────────────
ms <- load_built("MOLECULAR_SIGNATURES")
if (is.null(ms)) skip("database values valid", "MOLECULAR_SIGNATURES missing") else {
  assert_subset("database values valid", unique(as.character(cam$database)), names(ms))
  # PTMSIGDB is excluded by run_cameraPR()'s default: its sets carry a ";u"/";d" direction
  # suffix only the PTMsigDB branch of .prepare_sets() handles, and the released object has
  # no PTMSIGDB rows either.
  assert("PTMSIGDB not tested", !"PTMSIGDB" %in% as.character(cam$database))
}

# ── set metadata came from SET_TO_ID ──────────────────────────────────────────
s2i <- load_built("SET_TO_ID")
if (is.null(s2i)) skip("sets come from SET_TO_ID", "SET_TO_ID missing") else
  assert_subset("sets come from SET_TO_ID", unique(as.character(cam$set)),
                unique(as.character(s2i$set)))
for (cl in c("collection", "database", "set_id", "set"))
  assert(sprintf("CAMERA_RESULTS: %s fully joined", cl), !anyNA(cam[[cl]]),
         sprintf("%d NA", sum(is.na(cam[[cl]]))))

# ── contrast metadata came from CONTRAST_CONVERTER ────────────────────────────
cc <- load_built("CONTRAST_CONVERTER")
if (is.null(cc)) skip("contrasts come from CONTRAST_CONVERTER", "CONTRAST_CONVERTER missing") else {
  cc <- as.data.frame(cc)
  assert_subset("contrasts come from CONTRAST_CONVERTER",
                levels(droplevels(factor(cam$contrast))), as.character(cc$contrast))
  for (cl in c("contrast_type", "contrast_short"))
    assert(sprintf("CAMERA_RESULTS: %s fully joined", cl), !anyNA(cam[[cl]]),
           sprintf("%d NA", sum(is.na(cam[[cl]]))))
}

# ── statistics ────────────────────────────────────────────────────────────────
assert("CAMERA_RESULTS: p_value in [0,1]", all(cam$p_value >= 0 & cam$p_value <= 1, na.rm = TRUE))
assert("CAMERA_RESULTS: adj_p_value in [0,1]", all(cam$adj_p_value >= 0 & cam$adj_p_value <= 1, na.rm = TRUE))
bad <- sum(cam$adj_p_value < cam$p_value - 1e-9, na.rm = TRUE)
assert("CAMERA_RESULTS: adj_p_value >= p_value", bad == 0, sprintf("%d row(s) below", bad))
assert("CAMERA_RESULTS: direction is Up/Down",
       all(as.character(cam$direction) %in% c("Up", "Down")),
       paste(setdiff(unique(as.character(cam$direction)), c("Up", "Down")), collapse = ", "))
assert("CAMERA_RESULTS: t and z.std finite", all(is.finite(cam$t)) && all(is.finite(cam$z.std)),
       sprintf("%d non-finite t, %d non-finite z.std", sum(!is.finite(cam$t)), sum(!is.finite(cam$z.std))))

# set_size is the number of set members present in that dataset, so it cannot exceed the
# set's size in the database, and must clear run_cameraPR()'s min_size of 5
assert("CAMERA_RESULTS: set_size >= 5", all(cam$set_size >= 5L, na.rm = TRUE),
       sprintf("min = %s", min(cam$set_size, na.rm = TRUE)))
assert("CAMERA_RESULTS: set_size <= set_size_DB", all(cam$set_size <= cam$set_size_DB, na.rm = TRUE),
       sprintf("%d row(s) exceed", sum(cam$set_size > cam$set_size_DB, na.rm = TRUE)))
sr <- round(cam$set_size / cam$set_size_DB, 3L)
assert("CAMERA_RESULTS: size_ratio = set_size / set_size_DB",
       isTRUE(all.equal(sr, cam$size_ratio, tolerance = 1e-9)))

# BH is applied within tissue x assay x contrast x collection; re-derive one group and
# compare, which catches an adjustment done over the wrong grouping (or globally).
grp <- cam %>% dplyr::count(tissue, assay, contrast, collection) %>%
  dplyr::filter(n >= 20) %>% dplyr::slice(1)
if (nrow(grp) == 1) {
  g <- cam %>% dplyr::semi_join(grp, by = c("tissue", "assay", "contrast", "collection"))
  assert("CAMERA_RESULTS: BH adjusted within tissue x assay x contrast x collection",
         isTRUE(all.equal(p.adjust(g$p_value, method = "BH"), g$adj_p_value, tolerance = 1e-9)),
         sprintf("group %s/%s n=%d", grp$tissue[1], grp$assay[1], grp$n[1]))
} else skip("CAMERA_RESULTS: BH adjusted within tissue x assay x contrast x collection",
            "no group with >= 20 rows")

assert_unique("CAMERA_RESULTS: one row per tissue x assay x contrast x set x direction",
              paste(cam$tissue, cam$assay, cam$contrast, cam$set, cam$direction))

# ── diff vs the shipped package object (INFO only) ────────────────────────────
# Never a FAIL: this is computed from locally refit DA, so it moves with the freeze.
rda <- .find_rda("CAMERA_RESULTS")
if (is.na(rda)) report("diff vs package: CAMERA_RESULTS", "no package .rda") else {
  p <- as.data.frame(.load_rda(rda, "CAMERA_RESULTS"))
  report("diff vs package: CAMERA_RESULTS",
         sprintf("built %d x %d vs pkg %d x %d | cols %s | assays built [%s] pkg [%s]",
                 nrow(cam), ncol(cam), nrow(p), ncol(p),
                 if (identical(colnames(cam), colnames(p))) "identical" else "DIFFER",
                 paste(sort(unique(as.character(cam$assay))), collapse = ","),
                 paste(sort(unique(as.character(p$assay))), collapse = ",")))
}

finish()
