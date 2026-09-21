#!/usr/bin/env Rscript
# Step 03 — classify every freeze file against what the staging bucket already holds.
#
# Vendored from google_cloud_bucket_checks/03_diff_local_vs_staging.R. The adaptations
# are in the step README; the substantive ones are that the comparison is against the
# staging bucket rather than a production snapshot, and that a content change with no
# version bump is a CONFLICT instead of being silently skipped.
#
# Change types:
#   ADDED      freeze file whose stem is not in staging          -> upload
#   MODIFIED   different version, different md5                  -> upload, retire the old
#   REVERSIONED different version, same md5                      -> upload, retire the old; WARN
#   UNCHANGED  same version, same md5                            -> nothing to do
#   EQUIVALENT same version, different md5, same numbers         -> nothing to do; WARN
#   CONFLICT   same version, different md5                       -> FAIL
#   REPLACE    a CONFLICT, under ALLOW_CONFLICT=1                -> overwrite in place
#   DUPLICATE  staging holds two versions of one stem            -> reported, not touched
#   CARRIED    in staging, not regenerated, expected             -> reported, not touched
#   ORPHANED   in staging, not regenerated, not expected         -> reported, not touched

source(file.path(Sys.getenv("STAGE_LIB"), "bucket_helpers.R"))
source(file.path(Sys.getenv("STAGE_LIB"), "required_structure.R"))
source(file.path(.BH_ROOT, "scripts", "10_build_data", "lib", "file_versions.R"))

ALLOW_CONFLICT <- identical(Sys.getenv("ALLOW_CONFLICT"), "1")

# ---- Inputs ------------------------------------------------------------------------

snapshot_path <- file.path(SNAPSHOT_DIR, "latest_staging.tsv")
if (!file.exists(snapshot_path))
  stop("latest_staging.tsv not found. Run step 01 first: STEPS=01 make upload")

staging <- read.csv(snapshot_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                    colClasses = "character")
freeze  <- freeze_inventory()

message("freeze:  ", nrow(freeze), " upload artifacts under ", FREEZE_DIR)
message("staging: ", nrow(staging), " files at ", STAGING_BUCKET)

freeze$join_key  <- join_key(freeze$rel_path)
staging$join_key <- join_key(staging$rel_path)
freeze$version   <- vapply(freeze$rel_path,  strip_to_version, character(1), USE.NAMES = FALSE)
staging$version  <- vapply(staging$rel_path, strip_to_version, character(1), USE.NAMES = FALSE)

# ---- Resolve the incumbent in staging ----------------------------------------------
# A stem should appear once. When it appears more than once the bucket is carrying a
# stale leftover — a previous cycle uploaded a new version and failed to retire the old
# one. Compare against the highest version and report the rest.

staging <- staging[order(staging$join_key, -.ver_rank(staging$version)), , drop = FALSE]
is_incumbent <- !duplicated(staging$join_key)
duplicates   <- staging[!is_incumbent, , drop = FALSE]
incumbent    <- staging[is_incumbent, , drop = FALSE]

# ---- Classify each freeze file -----------------------------------------------------

rows <- lapply(seq_len(nrow(freeze)), function(i) {
  f   <- freeze[i, ]
  hit <- incumbent[incumbent$join_key == f$join_key, , drop = FALSE]
  meta <- parse_bic_name(f$rel_path)

  if (nrow(hit) == 0) {
    old_version <- NA_character_; old_md5 <- NA_character_; staging_gcs <- NA_character_
  } else {
    old_version <- hit$version[1]
    old_md5     <- hit$md5[1]
    staging_gcs <- hit$gcs_path[1]
  }

  # CONFLICT — same version string, different bytes. Upstream skips these outright, taking
  # an identical basename as proof of an identical file, but config/file_versions.json is
  # hand-maintained and this is exactly the mistake it can make: content regenerated,
  # version left alone. Uploading would overwrite a version already published under that
  # number; skipping would leave the bucket holding stale bytes under a number that now
  # means something else. Neither is safe to do silently, so the default is to stop and
  # make the caller choose. ALLOW_CONFLICT=1 chooses "overwrite in place" — the freeze is
  # right and the version is not moving. See classify_freeze_file() in lib/.
  change <- classify_freeze_file(f$version, old_version, f$md5, old_md5, ALLOW_CONFLICT)

  # Differing bytes are not yet a differing file. A recomputed numeric matrix can land on
  # the same values with different last bits, which md5 cannot see past, so before treating
  # this as a version mistake compare the numbers themselves. Only reached on CONFLICT and
  # REPLACE, so the download it costs happens on the rare path.
  if (change %in% c("CONFLICT", "REPLACE") &&
      freeze_numerically_equivalent(f$local_path, staging_gcs)) change <- "EQUIVALENT"

  data.frame(
    change_type      = change,
    rel_path         = f$rel_path,
    local_path       = f$local_path,
    gcs_path         = paste0(STAGING_BUCKET, "/", f$rel_path),
    staging_gcs_path = staging_gcs,
    join_key         = f$join_key,
    ome              = meta$ome,
    tissue_code      = meta$tissue_code,
    data_category    = meta$data_category,
    data_details     = meta$data_details,
    old_version      = old_version,
    new_version      = f$version,
    old_md5          = old_md5,
    new_md5          = f$md5,
    stringsAsFactors = FALSE
  )
})

diff_df <- do.call(rbind, rows)

# ---- Staging files the freeze does not regenerate ----------------------------------

orphan_keys <- setdiff(incumbent$join_key, freeze$join_key)
orphans <- incumbent[incumbent$join_key %in% orphan_keys, , drop = FALSE]

.staging_only_row <- function(df, change) {
  if (nrow(df) == 0) return(NULL)
  data.frame(
    change_type = change, rel_path = df$rel_path, local_path = NA_character_,
    gcs_path = NA_character_, staging_gcs_path = df$gcs_path, join_key = df$join_key,
    ome = NA_character_, tissue_code = NA_character_,
    data_category = NA_character_, data_details = NA_character_,
    old_version = df$version, new_version = NA_character_,
    old_md5 = df$md5, new_md5 = NA_character_,
    stringsAsFactors = FALSE
  )
}

diff_df <- rbind(
  diff_df,
  .staging_only_row(orphans[is_carried_forward(orphans$rel_path), , drop = FALSE], "CARRIED"),
  .staging_only_row(orphans[!is_carried_forward(orphans$rel_path), , drop = FALSE], "ORPHANED"),
  .staging_only_row(duplicates, "DUPLICATE")
)

# ---- Version drift against config/file_versions.json -------------------------------
# The freeze file names are stamped by freeze_version() at build time, so they agree with
# the map by construction. They stop agreeing when the map is edited without rebuilding,
# and then the version this step plans to upload is not the version the map decided on.
# A warning, not a failure: the fix is a rebuild, which is Stage 1's business.

drift <- vapply(seq_len(nrow(freeze)), function(i) {
  # Every freeze file is versioned from v2.0; a name with no token has no map entry to check.
  if (is.na(freeze$version[i])) return("")
  ext  <- sub("^\\.", "", path_ext(freeze$rel_path[i]))
  stem <- basename(join_key(freeze$rel_path[i]))
  mapped <- tryCatch(suppressMessages(freeze_version(stem, ext)), error = function(e) NA_character_)
  if (is.na(mapped) || identical(mapped, freeze$version[i])) "" else
    paste0(freeze$rel_path[i], " (built v", freeze$version[i], ", map says v", mapped, ")")
}, character(1))
drift <- drift[nzchar(drift)]

# ---- Write -------------------------------------------------------------------------

counts <- table(factor(diff_df$change_type,
                       levels = c("ADDED", "MODIFIED", "REVERSIONED", "REPLACE", "UNCHANGED",
                                  "EQUIVALENT", "CONFLICT", "DUPLICATE", "CARRIED",
                                  "ORPHANED")))
message("\nDiff summary:")
print(counts)

ts <- .timestamp()
diff_path <- file.path(DIFF_DIR, paste0("diff_", ts, ".tsv"))
report_write(diff_df, diff_path)
file.copy(diff_path, file.path(DIFF_DIR, "latest.tsv"), overwrite = TRUE)

status_for <- function(change, n) {
  if (n == 0) return("PASS")
  switch(change, CONFLICT = "FAIL", ORPHANED = "WARN", DUPLICATE = "WARN",
         REPLACE = "WARN", REVERSIONED = "WARN", EQUIVALENT = "WARN", "PASS")
}
report <- data.frame(
  status = vapply(names(counts), function(k) status_for(k, counts[[k]]), character(1)),
  check  = paste0("diff:", tolower(names(counts))),
  detail = paste0(as.integer(counts), " file(s)"),
  stringsAsFactors = FALSE
)
if (length(drift)) {
  report <- rbind(report, data.frame(
    status = "WARN", check = "diff:version-drift",
    detail = paste0(length(drift), " file(s) built at a version config/file_versions.json ",
                    "no longer maps to — rebuild Stage 1 or revert the map"),
    stringsAsFactors = FALSE
  ))
  message("\nVersion drift:"); message(paste(" ", drift, collapse = "\n"))
}
report_write(report, file.path(LOG_DIR, "upload_diff_report.tsv"))

# ---- Gate --------------------------------------------------------------------------

if (counts[["CONFLICT"]] > 0) {
  conflicts <- diff_df[diff_df$change_type == "CONFLICT", "rel_path"]
  message("\nCONFLICT — content changed but the version did not:")
  message(paste(" ", conflicts, collapse = "\n"))
  stop(length(conflicts), " file(s) differ from staging at the same version. ",
       "Bump them in config/file_versions.json and re-run Stage 1; or, if the freeze is ",
       "right and the version is not moving, re-run with ALLOW_CONFLICT=1 to overwrite ",
       "them in place. Both md5s are in ", diff_path)
}

if (counts[["REPLACE"]] > 0) {
  message("\nREPLACE — ALLOW_CONFLICT=1, so these will be OVERWRITTEN in place at their ",
          "current version rather than stopping the run:")
  message(paste(" ", diff_df[diff_df$change_type == "REPLACE", "rel_path"], collapse = "\n"))
}

if (counts[["REVERSIONED"]] > 0) {
  message("\nREVERSIONED — the version moved but the bytes did not. These upload and retire ",
          "the old object like any other version bump; the run does not stop. Check whether ",
          "the bump in config/file_versions.json was meant, since shipping one costs every ",
          "downstream consumer a version move for identical content:")
  message(paste(" ", diff_df[diff_df$change_type == "REVERSIONED", "rel_path"], collapse = "\n"))
}

if (counts[["ORPHANED"]] > 0) {
  message("\nORPHANED — in staging, not produced by this freeze, not on the carry-forward ",
          "list in lib/required_structure.R. Nothing is deleted; review them:")
  message(paste(" ", diff_df[diff_df$change_type == "ORPHANED", "rel_path"], collapse = "\n"))
}

message("Step 03 complete: ",
        counts[["ADDED"]] + counts[["MODIFIED"]] + counts[["REVERSIONED"]] + counts[["REPLACE"]],
        " file(s) to upload.")
