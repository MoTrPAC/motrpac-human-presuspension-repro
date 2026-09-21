#!/usr/bin/env Rscript
# Step 05 — validate the staging bucket against the manifest.
#
# Vendored from google_cloud_bucket_checks/05_validate_structure.R. Adaptations are in
# the step README; the two that change what this step does are that file matching is
# anchored to the composed file name rather than a set of unanchored greps, and that the
# value comparison against the installed data package is gone.
#
# Checks:
#   presence    every required manifest row resolves to a file in staging
#   columns     that file's header carries the required columns and none of the forbidden
#   forbidden   no version of any matched file carries a forbidden column, required row or not
#   coverage    every file in staging is claimed by a manifest row or the carry-forward list
#   qc-report   each subdir still has its HTML QC report
#
# Writes logs/structure_validation_<ts>.tsv and logs/upload_validate_report.tsv.
# Fails if any REQUIRED row fails.

source(file.path(Sys.getenv("STAGE_LIB"), "bucket_helpers.R"))
source(file.path(Sys.getenv("STAGE_LIB"), "required_structure.R"))
source(file.path(Sys.getenv("STAGE_LIB"), "expected_columns.R"))
# freeze_version() lives with Stage 1's map reader; step 03 reads it from the same place.
source(file.path(.BH_ROOT, "scripts", "10_build_data", "lib", "file_versions.R"))

FORCE_VALIDATE <- identical(Sys.getenv("FORCE_VALIDATE"), "1")

# ---- Is there anything to validate? ------------------------------------------------
# Step 04 is a dry run unless APPLY=1, so on a default run the bucket is exactly what it
# was before the stage started and validating it tells you about the PREVIOUS release,
# not this one. Skip rather than report failures that are really just "the upload has not
# happened yet". FORCE_VALIDATE=1 validates the bucket as it stands regardless, which is
# how you check a bucket someone else uploaded to.

diff_path <- file.path(DIFF_DIR, "latest.tsv")
if (!FORCE_VALIDATE && file.exists(diff_path)) {
  d <- read.csv(diff_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                colClasses = "character")
  # Every type step 04 uploads. REVERSIONED counts: its bytes are already in the bucket,
  # but at the old version, so the release is no more complete than for any other pending
  # upload until the new path exists.
  pending <- sum(d$change_type %in% c("ADDED", "MODIFIED", "REVERSIONED", "REPLACE"))
  if (pending > 0) {
    message("Step 03 still lists ", pending, " file(s) waiting to be uploaded, so the bucket ",
            "does not yet hold this release.\nRun step 04 with APPLY=1 first, or set ",
            "FORCE_VALIDATE=1 to validate the bucket as it stands.")
    quit(save = "no", status = 77)
  }
}

# Listed fresh rather than read from step 01's snapshot: step 04 changes the bucket, and a
# snapshot taken before it ran would have this step validating files that are no longer
# there and missing the ones that are.
staging <- gcs_inventory(STAGING_BUCKET)
staging$rel_path <- sub(paste0("^", STAGING_BUCKET, "/"), "", staging$gcs_path)
message("Validating ", nrow(staging), " files at ", STAGING_BUCKET)

# ---- Matching ----------------------------------------------------------------------
# Upstream greps the whole listing once per manifest field and keeps whatever survives all
# five. Unanchored substring matching is loose enough to cross rows — "da" appears inside
# "metadata", "samples" inside "removed-samples" (which is why upstream has to wrap
# data_details in underscores) — so this composes the file name the manifest describes and
# matches it as a prefix. A file that does not sit where its own name says it should is a
# miss here, and coverage below reports it.

.category_dir <- function(data_category) {
  if (data_category %in% c("qc-norm", "imputed")) "qc-norm" else data_category
}

.expected_prefix <- function(row) {
  paste0(row$gcs_subdir, "/", .category_dir(row$data_category), "/",
         FILE_HEADER, "_", row$tissue_code, "_", row$ome, "_",
         row$data_category, "_", row$data_details, "_v")
}

# First line of a GCS object, read as a range request. Upstream downloads whole files;
# a qc-norm matrix is hundreds of MB and only its header is under test here.
#
# Column names are unquoted, except in the methylcap DA tables — those are external MALAX
# GLMM output copied into the freeze verbatim, and they were written with quote = TRUE.
# read.csv() would strip the quotes for free; splitting the raw line does not, so it is
# done here. Without it every methylcap column reads as "feature_id" with the quotes
# attached and the whole schema looks absent.
.read_header <- function(gcs_path) {
  out <- suppressWarnings(system(paste(GSUTIL, "cat -r 0-1048575", shQuote(gcs_path)),
                                 intern = TRUE, ignore.stderr = TRUE))
  if (!length(out)) return(character(0))
  gsub('^"|"$', "", strsplit(out[1], "\t", fixed = TRUE)[[1]])
}

# Every file is versioned from v2.0, but a pre-v2.0 incumbent in the bucket may carry no
# token, and "found at vNA" would read as a parse failure rather than a file with none.
.found_at <- function(rel_path) {
  v <- strip_to_version(rel_path)
  if (is.na(v)) "file found (unversioned)" else paste0("file found at v", v)
}

.result <- function(row, status, reason, file_path = NA_character_) {
  data.frame(
    ome = row$ome, tissue = row$tissue, tissue_code = row$tissue_code,
    data_category = row$data_category, data_details = row$data_details,
    required = row$required, status = status, reason = reason,
    file_path = file_path, stringsAsFactors = FALSE
  )
}

matched_paths <- character(0)

# Files in the bucket carrying a version config/file_versions.json does not name for them.
# Step 03 checks the same map against the FREEZE, before anything is uploaded; this checks
# it against what the bucket actually ended up holding, which is the direction a partial or
# reordered upload breaks. Collected here and reported once, rather than folded into the
# per-row result, because a version mismatch and a column mismatch are different faults.
version_mismatches <- character(0)

# Files in the bucket whose header carries a column the schema forbids. Checked on every
# version of a stem present, not only the one that would be read, and whether or not the
# schema requires anything: a removed-samples file requires no column but must never carry
# a participant identifier. Case-insensitive. Any hit fails the run, required row or optional.
forbidden_hits <- character(0)

.check_forbidden <- function(rel_paths, forbidden, headers = list()) {
  for (p in rel_paths) {
    if (grepl("\\.gz$", p)) next
    header <- if (!is.null(headers[[p]])) headers[[p]] else .read_header(paste0(STAGING_BUCKET, "/", p))
    extra  <- forbidden[tolower(forbidden) %in% tolower(header)]
    if (length(extra))
      forbidden_hits <<- c(forbidden_hits, paste0(p, " (", paste(extra, collapse = ", "), ")"))
  }
  invisible(NULL)
}

.check_version <- function(rel_path) {
  v <- strip_to_version(rel_path)
  if (is.na(v)) return(invisible(NULL))   # a pre-v2.0 incumbent with no version token
  ext    <- sub("^\\.", "", path_ext(rel_path))
  stem   <- basename(join_key(rel_path))
  mapped <- tryCatch(suppressMessages(freeze_version(stem, ext)), error = function(e) NA_character_)
  if (is.na(mapped) || identical(mapped, v)) return(invisible(NULL))
  version_mismatches <<- c(version_mismatches,
                           paste0(rel_path, " (bucket v", v, ", map says v", mapped, ")"))
  invisible(NULL)
}

validate_row <- function(row) {
  prefix  <- .expected_prefix(row)
  matched <- staging$rel_path[startsWith(staging$rel_path, prefix)]

  if (length(matched) == 0) {
    if (!row$required)
      return(.result(row, "PASS", paste0("no file — optional (", row$data_details, ")")))
    return(.result(row, "FAIL", paste0("no file matching ", prefix, "*")))
  }

  # Highest version wins when a stem is present more than once. Step 03 reports the
  # duplicate; this step validates the one that would be read.
  vers <- vapply(matched, strip_to_version, character(1), USE.NAMES = FALSE)
  best <- matched[which.max(.ver_rank(vers))]
  matched_paths <<- c(matched_paths, matched)
  .check_version(best)

  expected <- .get_expected_cols(row$data_category, row$ome, row$data_details)
  if (is.null(expected))
    return(.result(row, "PASS", paste0(.found_at(best), "; no column schema defined"), best))

  # Superseded versions still in the bucket: forbidden columns only, since they are not read.
  if (length(expected$forbidden))
    .check_forbidden(setdiff(matched, best), expected$forbidden)

  # A range read returns gzip bytes, not text. Nothing ships compressed today, so rather
  # than teach the reader to decompress, say plainly that the schema went unchecked —
  # "header could not be read" would look like a broken file.
  if (grepl("\\.gz$", best))
    return(.result(row, "PASS", paste0(.found_at(best), "; compressed, schema not checked"), best))

  header <- .read_header(paste0(STAGING_BUCKET, "/", best))
  if (!length(header))
    return(.result(row, "FAIL", "file found in staging but its header could not be read", best))
  if (length(expected$forbidden))
    .check_forbidden(best, expected$forbidden, stats::setNames(list(header), best))
  if (length(expected$required) == 0 && length(expected$forbidden) == 0)
    return(.result(row, "PASS", paste0(.found_at(best), "; no column schema defined"), best))

  missing <- setdiff(expected$required, header)
  extra   <- intersect(expected$forbidden, header)

  # A qc-norm or imputed matrix is feature_id plus one column per sample; a header of
  # feature_id alone is a matrix that lost its samples.
  if (row$data_category %in% c("qc-norm", "imputed") && length(header) < 2)
    return(.result(row, "FAIL", "matrix has no sample columns", best))

  if (length(missing) || length(extra)) {
    parts <- c(if (length(missing)) paste("missing cols:", paste(missing, collapse = ", ")),
               if (length(extra))   paste("forbidden cols present:", paste(extra, collapse = ", ")))
    return(.result(row, "FAIL", paste(parts, collapse = "; "), best))
  }

  .result(row, "PASS", paste0(.found_at(best), "; ", length(header), " columns"), best)
}

results <- do.call(rbind, lapply(seq_len(nrow(required_structure)), function(i) {
  validate_row(required_structure[i, ])
}))

# ---- resources/ --------------------------------------------------------------------

resource_results <- do.call(rbind, lapply(seq_len(nrow(REQUIRED_RESOURCES)), function(i) {
  r      <- REQUIRED_RESOURCES[i, ]
  prefix <- paste0(r$rel_stem, if (r$versioned) "_v" else ".")
  hit    <- staging$rel_path[startsWith(staging$rel_path, prefix)]
  matched_paths <<- c(matched_paths, hit)
  row <- data.frame(ome = NA_character_, tissue = NA_character_, tissue_code = NA_character_,
                    data_category = "resources", data_details = basename(r$rel_stem),
                    required = r$required, stringsAsFactors = FALSE)
  if (length(hit) == 0) .result(row, if (r$required) "FAIL" else "PASS",
                                paste0("no file matching ", prefix, "*"))
  else .result(row, "PASS", .found_at(hit[1]), hit[1])
}))

results <- rbind(results, resource_results)

# ---- Coverage ----------------------------------------------------------------------
# Every file in the bucket should be claimed by a manifest row or be a known carry-forward.
# Anything left over is a file nothing in this pipeline knows about: a stale artifact from
# a retired layout, or a new one that belongs in lib/required_structure.R.

unclaimed <- setdiff(staging$rel_path, unique(matched_paths))
unclaimed <- unclaimed[!is_carried_forward(unclaimed)]

# ---- QC reports --------------------------------------------------------------------

missing_qc <- QC_REPORT_SUBDIRS[!vapply(QC_REPORT_SUBDIRS, function(s)
  any(grepl(paste0("^", s, "/.*qc-report.*\\.html$"), staging$rel_path)), logical(1))]

# ---- Write -------------------------------------------------------------------------

ts <- .timestamp()
val_log <- file.path(LOG_DIR, paste0("structure_validation_", ts, ".tsv"))
report_write(results, val_log)

message("\nValidation summary:")
print(table(results$status))

required_fails <- results[results$status == "FAIL" & results$required, , drop = FALSE]
optional_fails <- results[results$status == "FAIL" & !results$required, , drop = FALSE]

report <- data.frame(
  status = c(if (nrow(required_fails)) "FAIL" else "PASS",
             if (nrow(optional_fails)) "WARN" else "PASS",
             if (length(unclaimed))    "WARN" else "PASS",
             if (length(missing_qc))   "WARN" else "PASS",
             if (length(version_mismatches)) "FAIL" else "PASS",
             if (length(forbidden_hits)) "FAIL" else "PASS"),
  check  = c("validate:required", "validate:optional", "validate:coverage", "validate:qc-report",
             "validate:versions", "validate:forbidden"),
  detail = c(paste0(sum(results$required) - nrow(required_fails), "/", sum(results$required),
                    " required entries present and well-formed"),
             paste0(nrow(optional_fails), " optional entries failed"),
             paste0(length(unclaimed), " staging file(s) claimed by no manifest row"),
             if (length(missing_qc)) paste("no HTML QC report in:", paste(missing_qc, collapse = ", "))
             else "every subdir has an HTML QC report",
             if (length(version_mismatches))
               paste0(length(version_mismatches), " file(s) at a version config/file_versions.json ",
                      "does not name for them")
             else "every versioned file carries the version the map decided",
             if (length(forbidden_hits))
               paste0(length(forbidden_hits), " file(s) carry a forbidden column")
             else "no file carries a forbidden column"),
  stringsAsFactors = FALSE
)
report_write(report, file.path(LOG_DIR, "upload_validate_report.tsv"))

if (length(unclaimed)) {
  message("\nStaging files claimed by no manifest row:")
  message(paste(" ", unclaimed, collapse = "\n"))
}
if (length(missing_qc))
  message("\nWARN: no HTML QC report found in: ", paste(missing_qc, collapse = ", "))
if (length(forbidden_hits)) {
  message("\nFORBIDDEN COLUMNS present:")
  message(paste(" ", forbidden_hits, collapse = "\n"))
}
if (length(version_mismatches)) {
  message("\nFILE VERSIONS the map does not name:")
  message(paste(" ", version_mismatches, collapse = "\n"))
}

if (nrow(required_fails) > 0) {
  message("\nFailed required entries:")
  print(required_fails[, c("ome", "tissue_code", "data_category", "data_details", "reason")])
  stop(nrow(required_fails), " required manifest entries FAILED. Inspect ", val_log)
}

# The driver reads this step's exit code, not its report, so a FAIL row has to stop here to
# mean anything.
if (length(version_mismatches)) {
  stop(length(version_mismatches), " bucket file(s) carry a version config/file_versions.json ",
       "does not name for them. Either the upload put the wrong version in the bucket, or the ",
       "map moved after the upload; the listing above names both versions for each.")
}

if (length(forbidden_hits)) {
  stop(length(forbidden_hits), " bucket file(s) carry a column their schema forbids; the ",
       "listing above names each file and column. Remove or replace them before promoting.")
}

message("Step 05 complete.")
