#!/usr/bin/env Rscript
# Tests for the qc-norm FREEZE (step 06). Downstream: step 07 (feature map) and step 08
# (the *_QC objects) read these txt files; each must be complete, single-version, and
# have the right columns. Diff each qc-norm matrix against the reference bucket.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
suppressWarnings(suppressMessages(library(data.table)))
cat("== qc-norm freeze tests ==\n")

# 1. no stale duplicate-version freeze files (the v1.2/v1.4 bug this suite exists to catch)
assert_no_dup_versions()

freeze <- list.files(.T_FREEZE, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
assert("freeze files present", length(freeze) > 0, "staging/freeze empty — run step 06 first")
if (length(freeze) == 0) finish()

# 1b. Manifest presence: EVERY row of OME_TISSUE_CODE must have its generated files, and the
#     qc-norm file must carry the data_details string that ome actually uses. Ported from
#     data-raw/google_cloud_bucket_checks/required_structure.R, which builds the same manifest
#     (its `ome_tissue_pairs`) and drives the bucket validator's file-presence check.
#
#     One deliberate difference from upstream's filter: the `lab-*` omes are excluded here,
#     because this pipeline's OME_TISSUE_CODE is built from the quant-id catalog and carries
#     5 clinical-chemistry lab rows that upstream's copy does not have — they have no
#     qc_details mapping and produce no freeze files. That leaves 48 rows.
#
#     Checking data_details matters: it is what distinguishes a correctly named matrix from
#     one written with the wrong transform (log2 vs log-cpm vs beta-values vs absolute).
#     Upstream's own manifest turns an unmapped ome into the literal search term "_NA_", so the
#     mapping is asserted here rather than allowed to silently miss.
otc <- load_preflight("OME_TISSUE_CODE")
if (!is.null(otc)) {
  qc_details_for <- function(ome) {
    if (ome == "transcript-rna-seq")                  return("log-cpm")
    if (ome %in% c("prot-pr", "prot-ph"))             return("log2-mn")
    if (ome == "prot-ol")                             return("log2")
    if (ome == "epigen-atac-seq")                     return("log-cpm")
    if (ome == "epigen-methylcap-seq")                return("beta-values")
    # must precede the metab catch-all, as upstream notes
    if (ome %in% c("prot-clinical", "metab-t-clinical")) return("absolute")
    if (grepl("^metab", ome))                         return("log2")
    return(NA_character_)
  }
  man <- otc[!grepl("^lab-", otc$ome), ]
  man$token <- paste0(man$tissue_code, "_", man$ome)
  man$qc_details <- vapply(man$ome, qc_details_for, character(1), USE.NAMES = FALSE)

  assert("every manifest ome maps to a qc-norm data_details",
         all(!is.na(man$qc_details)),
         sprintf("unmapped: %s", paste(unique(man$ome[is.na(man$qc_details)]), collapse = ", ")))

  fz_base <- basename(freeze)
  missing_for <- function(patterns)
    man$token[!vapply(patterns, function(p) any(grepl(p, fz_base, fixed = TRUE)), logical(1),
                      USE.NAMES = FALSE)]

  check_present <- function(desc, patterns) {
    miss <- missing_for(patterns)
    assert(desc, length(miss) == 0,
           sprintf("missing %d of %d: %s", length(miss), nrow(man),
                   paste(utils::head(miss, 5), collapse = ", ")))
  }
  check_present("every OME_TISSUE_CODE row has a qc-norm file",
                paste0(man$token, "_qc-norm_", man$qc_details))
  check_present("every OME_TISSUE_CODE row has metadata_features",
                paste0(man$token, "_metadata_features"))
  check_present("every OME_TISSUE_CODE row has metadata_samples",
                paste0(man$token, "_metadata_samples"))

  # removed-samples is required = FALSE upstream: an ome x tissue with no dropped samples
  # legitimately has no file, so this is reported rather than asserted.
  rm_miss <- missing_for(paste0(man$token, "_metadata_removed-samples"))
  report("removed-samples coverage",
         sprintf("%d of %d manifest rows have a file", nrow(man) - length(rm_miss), nrow(man)))
} else skip("freeze manifest presence", "OME_TISSUE_CODE preflight object missing")

tok_of <- function(p, type) sub(paste0("_", type, ".*"), "",
                                sub("^human-precovid-sed-adu_", "", basename(p)))
ome_of <- function(p, type) sub("^[^_]+_", "", tok_of(p, type))

# 2. Column schema for EVERY file, per ome. Ported from
#    data-raw/google_cloud_bucket_checks/expected_columns.R, which the bucket validator
#    applies to every manifest row (05_validate_structure.R:130-152). This previously checked
#    three generic things on head(files, 3), so 45 of 48 files of each kind were never
#    schema-checked at all and a metadata_features file missing gene_symbol / uniprot /
#    refmet_name / ensembl_gene / entrez_gene / relationship_to_gene passed silently.
#
#    Upstream's `forbidden` column mechanism is not ported: it is empty for all 22 of its
#    schemas. The identifier check on removed-samples (2b) is this repo's own.
qc_files <- freeze[grepl("qc-norm", basename(freeze))]

# metadata_features requirements resolve by ome prefix, matching upstream's fallback order.
# Note prot-clinical lands on the prot set and metab-t-clinical on the metab set, which is
# what upstream does too (expected_columns.R:86-88).
features_required_for <- function(ome) {
  if (ome == "transcript-rna-seq") return(c("assay", "feature_id", "gene_symbol", "ensembl_gene"))
  if (grepl("^prot", ome))         return(c("assay", "feature_id", "gene_symbol", "uniprot"))
  if (grepl("^metab", ome))        return(c("assay", "feature_id", "refmet_name"))
  if (grepl("^epigen", ome))       return(c("assay", "feature_id", "gene_symbol", "ensembl_gene",
                                            "entrez_gene", "relationship_to_gene"))
  return(character(0))
}
# fread, not read.csv, in this file only: these checks need a header, one column, or 20
# rows from each matrix, and read.csv cannot read part of a file.
header_of <- function(f) names(fread(f, nrows = 0, sep = "\t", header = TRUE))

for (f in qc_files) {
  # Type the value columns numeric instead of letting fread infer them. A sample column with
  # no measurement at all is every-value-NA, and a bare NA infers as logical — not a type
  # error, just an empty sample. prot-clinical has one (vial 11225060213, NA for
  # CK/Glucagon/Insulin) and the reference release carries it too. Forcing numeric keeps those
  # NAs numeric; fread declines to downgrade a genuine string column, so the check below still
  # catches real corruption.
  n_columns <- length(header_of(f))
  d <- fread(f, nrows = 20, sep = "\t", header = TRUE, data.table = FALSE,
             colClasses = c("character", rep("numeric", n_columns - 1)))
  b <- basename(f)
  assert(sprintf("qc-norm %s: has feature_id", b), names(d)[1] == "feature_id",
         sprintf("first column is %s", names(d)[1]))
  # upstream's expected_columns.R:26 says sample columns are "checked by nrow/ncol"; no such
  # check exists anywhere in that folder, so a matrix with no samples at all passes it.
  assert(sprintf("qc-norm %s: has >= 1 sample column", b), ncol(d) > 1,
         sprintf("only %d column(s)", ncol(d)))
  is_numeric_column <- vapply(d[, -1, drop = FALSE], is.numeric, logical(1))
  assert(sprintf("qc-norm %s: value columns numeric", b),
         ncol(d) > 1 && all(is_numeric_column),
         sprintf("non-numeric: %s",
                 paste(names(is_numeric_column)[!is_numeric_column], collapse = ", ")))
}
for (f in freeze[grepl("metadata_features", basename(freeze))]) {
  req <- features_required_for(ome_of(f, "metadata_features"))
  if (length(req) == 0) { skip(sprintf("features %s: required columns", basename(f)), "no schema for this ome"); next }
  h <- header_of(f)
  assert(sprintf("features %s: required columns", basename(f)), all(req %in% h),
         sprintf("missing: %s", paste(setdiff(req, h), collapse = ", ")))
}
for (f in freeze[grepl("metadata_samples", basename(freeze))]) {
  h <- header_of(f)
  assert(sprintf("samples %s: has vialLabel", basename(f)), "vialLabel" %in% h)
}

# 2b. removed-samples files carry no participant identifiers, demographics or study design.
#     Checked on the vendored sources as well as the staged copies, so a re-fetch from the
#     bucket fails here even when STEMS skips stage_removed_samples.
identifying_columns <- c("pid", "sex", "age", "race", "siteid", "short.id",
                         "rdgroup", "visit", "timept")
removed_sources <- list.files(
  file.path(.T_ROOT, "scripts", "00_preflight", "data-raw", "sources", "removed_samples"),
  pattern = "removed-samples.*\\.txt$", full.names = TRUE)
assert("removed-samples sources present", length(removed_sources) > 0)
for (f in c(removed_sources, freeze[grepl("metadata_removed-samples", basename(freeze))])) {
  present <- intersect(tolower(header_of(f)), identifying_columns)
  assert(sprintf("removed-samples %s: no identifying columns (%s)",
                 basename(f), basename(dirname(dirname(f)))),
         length(present) == 0,
         sprintf("present: %s", paste(present, collapse = ", ")))
}

# 3. the ids in each metadata file must correspond exactly to the axes of its qc-norm
#    matrix: metadata_samples$vialLabel to the value columns, metadata_features$feature_id
#    to the feature column. The checks above only confirm that files exist and carry the
#    right headers, so a matrix written with a stale sample list — or one left truncated
#    by an interrupted write — passes them and is caught here. Uniqueness is asserted in the
#    same pass, on the vectors already read, since these matrices run to hundreds of MB.
path_by_token <- function(type) {
  p <- freeze[grepl(type, basename(freeze), fixed = TRUE)]
  stats::setNames(p, vapply(p, tok_of, "", type = type))
}
samp_paths <- path_by_token("metadata_samples")
feat_paths <- path_by_token("metadata_features")

# one assertion per axis, reporting both directions: ids documented but missing from the
# matrix, and ids in the matrix with no metadata row.
assert_ids_match <- function(desc, meta_ids, axis_ids, meta_label, axis_label) {
  absent <- setdiff(meta_ids, axis_ids)
  undocumented <- setdiff(axis_ids, meta_ids)
  assert(desc, length(absent) == 0 && length(undocumented) == 0,
         sprintf("%d %s absent from %s (e.g. %s); %d %s undocumented (e.g. %s)",
                 length(absent), meta_label, axis_label,
                 paste(utils::head(absent, 3), collapse = ", "),
                 length(undocumented), axis_label,
                 paste(utils::head(undocumented, 3), collapse = ", ")))
}
# read only what each check needs — these matrices run to hundreds of MB
first_col <- function(f) as.character(fread(f, select = 1, sep = "\t", header = TRUE)[[1]])
id_col <- function(f, col) as.character(fread(f, select = col, sep = "\t", header = TRUE)[[1]])

for (f in qc_files) {
  tk <- tok_of(f, "qc-norm")
  qc_samples <- names(fread(f, nrows = 0, sep = "\t", header = TRUE))[-1]
  qc_features <- first_col(f)

  # A duplicate id is never legitimate here and silently fans out any downstream join.
  # Neither this suite nor the bucket validator checked for it.
  assert_unique(sprintf("%s: qc-norm feature_ids unique", tk), qc_features)
  assert_unique(sprintf("%s: qc-norm sample columns unique", tk), qc_samples)

  d <- sprintf("%s: metadata_samples ↔ qc-norm columns", tk)
  if (tk %in% names(samp_paths)) {
    sm_ids <- id_col(samp_paths[[tk]], "vialLabel")
    assert_unique(sprintf("%s: metadata_samples vialLabels unique", tk), sm_ids)
    assert_ids_match(d, sm_ids, qc_samples, "metadata sample(s)", "qc-norm columns")
  } else skip(d, "no metadata_samples file for this token")

  d <- sprintf("%s: metadata_features ↔ qc-norm rows", tk)
  if (tk %in% names(feat_paths)) {
    fm_ids <- id_col(feat_paths[[tk]], "feature_id")
    assert_unique(sprintf("%s: metadata_features feature_ids unique", tk), fm_ids)
    assert_ids_match(d, fm_ids, qc_features, "metadata feature(s)", "qc-norm rows")
  } else skip(d, "no metadata_features file for this token")
}

# 4. diff each qc-norm matrix vs the reference bucket (SKIPs cleanly without gsutil)
for (f in qc_files) diff_freeze_vs_bucket(f)
finish()
