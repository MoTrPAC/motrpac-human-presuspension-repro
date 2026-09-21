#!/usr/bin/env Rscript
# Tests for the *_QC objects (step 08). Downstream: load_qc / the DA generators index
# each leaf as list(qc_norm, sample_metadata, feature_metadata); qc_norm must be a
# numeric data.frame with feature rownames and vialLabel colnames; sample_metadata must
# carry the design + per-ome covariate columns; blood transcriptomics is split 450/rest.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
cat("== qc object tests ==\n")

# 1. assemble the nested [[tissue]][[ome]] list with load_qc_local() — the loader the DA
#    generators use — called with ITS OWN defaults, so the tests see exactly the filtering a
#    downstream consumer gets: load_acute_only = TRUE, remove_redundant_metab = TRUE,
#    remove_unnamed_metab = FALSE, load_clinical = FALSE.
#    load_qc_local lives in lib/qc_helpers.R alongside the other QC-object helpers.
source(file.path(.root, "scripts", "10_build_data", "lib", "qc_helpers.R"))

qc <- tryCatch(load_qc_local(), error = function(e) structure(NULL, msg = conditionMessage(e)))
assert("load_qc_local assembles the nested list", !is.null(qc), attr(qc, "msg") %||% "returned NULL")

if (!is.null(qc)) {
  # 2. completeness — every expected tissue×ome present. include_clinical tracks
  #    load_qc_local's load_clinical default: the clinical omes are filtered out by default,
  #    so counting them here would read as missing. Section 3b covers them instead.
  present <- unlist(lapply(names(qc), function(t) paste0(t, ".", names(qc[[t]]))))
  assert_completeness("all expected tissue×ome present", present, expected_tissue_omes(include_clinical = FALSE))

  # 3. per-leaf contract (structure, classes, alignment, required + covariate columns)
  for (t in names(qc)) for (o in names(qc[[t]])) assert_qc_leaf(t, o, qc[[t]][[o]])
}

# 3b. step 08 builds the clinical omes, but load_qc_local's defaults filter them out, so they
#     need an explicit opt-in pass to stay covered. selected_omes keeps this to the two
#     objects rather than re-loading all 47.
clin_omes <- c("prot-clinical", "metab-t-clinical")
qc_clin <- tryCatch(load_qc_local(selected_omes = clin_omes, load_clinical = TRUE),
                    error = function(e) structure(NULL, msg = conditionMessage(e)))
if (!is.null(qc_clin)) {
  for (t in names(qc_clin)) for (o in intersect(names(qc_clin[[t]]), clin_omes))
    assert_qc_leaf(t, o, qc_clin[[t]][[o]])
} else skip("clinical leaf contract", "load_qc_local(load_clinical = TRUE) failed")

# 4. blood transcriptomics split invariants (raw objects, before recombination)
b1 <- load_built("blood_transcript_1"); b2 <- load_built("blood_transcript_2")
if (!is.null(b1) && !is.null(b2) && !is.null(b1$qc_norm) && !is.null(b2$qc_norm)) {
  assert("blood_transcript_1 qc_norm has 450 columns", ncol(b1$qc_norm) == 450, sprintf("ncol=%d", ncol(b1$qc_norm)))
  assert("blood split feature rownames identical (cbind-safe)", identical(rownames(b1$qc_norm), rownames(b2$qc_norm)))
  assert("blood split samples disjoint", length(intersect(colnames(b1$qc_norm), colnames(b2$qc_norm))) == 0)

  # Sample coverage. The three checks above all still pass if the tail of the matrix is
  # dropped outright — 450 columns, matching rownames and disjointness say nothing about
  # what happened to columns 451:n. These pin the split to the sample list instead.
  # qc_norm_results.R copies the whole BLOOD_TRNSCRPT leaf into both halves and replaces
  # only $qc_norm, so either half's sample_metadata is the full post-merge (pheno-joined)
  # table and is the right thing to measure coverage against.
  if (!is.null(b1$sample_metadata) && !is.null(b1$sample_metadata$vialLabel)) {
    meta_vials  <- unique(as.character(b1$sample_metadata$vialLabel))
    split_vials <- union(colnames(b1$qc_norm), colnames(b2$qc_norm))

    assert_completeness("blood split covers every merged-metadata sample", split_vials, meta_vials)

    extra <- setdiff(split_vials, meta_vials)
    assert("blood split has no sample absent from merged metadata", length(extra) == 0,
           sprintf("%d extra: %s", length(extra), paste(utils::head(extra, 8), collapse = ", ")))

    # Arithmetic cross-check on the union: it would catch an overlap between the halves
    # independently of the disjointness assert above.
    assert("blood split column counts sum to the covered sample count",
           ncol(b1$qc_norm) + ncol(b2$qc_norm) == length(split_vials),
           sprintf("%d + %d != %d", ncol(b1$qc_norm), ncol(b2$qc_norm), length(split_vials)))

    assert("blood split halves carry the same sample_metadata",
           identical(b1$sample_metadata, b2$sample_metadata),
           "halves disagree — the coverage check above measures against half 1 only")
  } else skip("blood split sample coverage", "sample_metadata or its vialLabel column missing")
} else skip("blood split invariants", "blood_transcript_1/2 not built or missing qc_norm")

# 5. diff each *_QC vs the Data package (feature/sample overlap + per-feature correlation)
# diff each *_QC's qc_norm vs the package (REPORT only — a reproduction difference is
# logged, not a build-stopping FAIL; assert_qc_leaf above owns the hard structural checks).
report_qc_cor <- function(d, built_qc, pkg_qc) {
  if (is.null(built_qc) || is.null(pkg_qc)) return(report(d, "missing qc_norm"))
  cf <- intersect(rownames(built_qc), rownames(pkg_qc))
  cs <- intersect(colnames(built_qc), colnames(pkg_qc))
  if (length(cf) < 3 || length(cs) < 3) return(report(d, "few shared features/samples"))
  med <- stats::median(.per_feat_cor(as.matrix(built_qc[cf, cs]),
                                     as.matrix(pkg_qc[cf, cs])), na.rm = TRUE)
  report(sprintf("%s median per-feat cor", d), sprintf("%.4f", med))
}
diff_one_qc <- function(nm) {
  d <- sprintf("diff qc vs package: %s", nm)
  rda <- .find_rda(nm); if (is.na(rda)) return(report(d, "no package .rda"))
  report_qc_cor(d, load_built(nm)$qc_norm, .load_rda(rda, nm)$qc_norm)
}
for (nm in sub("\\.rda$", "", list.files(.T_DATA, pattern = "_QC\\.rda$"))) diff_one_qc(nm)

# Blood transcriptomics is not named *_QC — it ships split 450/rest, here and in the package
# alike — so the loop above never reached it, leaving it the one ome with no object-level
# diff. Recombine both sides the way load_qc() does (cbind half 2 onto half 1) and compare
# that. cbind pairs data.frame rows by position, so mismatched rownames would silently
# mis-pair features; the halves are checked for that above on the built side, and below on
# the package side.
local({
  description   <- "diff qc vs package: blood_transcript_1+2"
  half_names    <- c("blood_transcript_1", "blood_transcript_2")
  package_files <- vapply(half_names, .find_rda, character(1))
  if (anyNA(package_files)) return(report(description, "no package .rda"))

  built_halves   <- lapply(half_names, load_built)
  package_halves <- Map(.load_rda, package_files, half_names)
  qc_norm_of <- function(half) if (is.null(half)) NULL else half$qc_norm

  if (any(vapply(c(built_halves, package_halves),
                 function(half) is.null(qc_norm_of(half)), logical(1))))
    return(report(description, "missing qc_norm"))
  if (!identical(rownames(qc_norm_of(package_halves[[1]])),
                 rownames(qc_norm_of(package_halves[[2]]))))
    return(report(description, "package halves have mismatched feature rownames — cbind unsafe"))

  report_qc_cor(description,
                cbind(qc_norm_of(built_halves[[1]]),   qc_norm_of(built_halves[[2]])),
                cbind(qc_norm_of(package_halves[[1]]), qc_norm_of(package_halves[[2]])))
})
finish()
