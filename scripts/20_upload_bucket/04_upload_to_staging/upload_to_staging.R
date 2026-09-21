#!/usr/bin/env Rscript
# Step 04 — apply the step-03 diff to the staging bucket.
#
# Vendored from google_cloud_bucket_checks/04_upload_to_staging.R. Same shape: upload the
# ADDED, MODIFIED and REVERSIONED rows, verify each upload's md5 against the local file, and
# only then retire the version it supersedes. Adaptations are in the step README.
#
# Writes logs/upload_log_<ts>.tsv and logs/upload_apply_report.tsv.

source(file.path(Sys.getenv("STAGE_LIB"), "bucket_helpers.R"))

APPLY <- identical(Sys.getenv("APPLY"), "1")

diff_path <- file.path(DIFF_DIR, "latest.tsv")
if (!file.exists(diff_path))
  stop("diffs/latest.tsv not found. Run step 03 first: STEPS=03 make upload")

diff_df  <- read.csv(diff_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                     colClasses = "character")
to_upload <- diff_df[diff_df$change_type %in% c("ADDED", "MODIFIED", "REVERSIONED", "REPLACE"), ,
                     drop = FALSE]

# The set that supersedes an object at another version, and so retires it after verifying.
# REVERSIONED is here for the same reason MODIFIED is: the version moved, so the upload
# lands at a new path and the old one would otherwise linger as a DUPLICATE next cycle.
SUPERSEDING <- c("MODIFIED", "REVERSIONED")

message("Files to upload: ", nrow(to_upload))
message("  ADDED:       ", sum(to_upload$change_type == "ADDED"))
message("  MODIFIED:    ", sum(to_upload$change_type == "MODIFIED"),
        " (the superseded version is removed from staging after verification)")
message("  REVERSIONED: ", sum(to_upload$change_type == "REVERSIONED"),
        " (same bytes at a new version — uploaded and superseded like MODIFIED; step 03 warned)")
message("  REPLACE:     ", sum(to_upload$change_type == "REPLACE"),
        " (overwritten in place at the same version — step 03 ran with ALLOW_CONFLICT=1)")

if (nrow(to_upload) == 0) {
  report_write(data.frame(status = "PASS", check = "upload:nothing-to-do",
                          detail = "staging already matches the freeze",
                          stringsAsFactors = FALSE),
               file.path(LOG_DIR, "upload_apply_report.tsv"))
  message("Step 04 complete: nothing to upload.")
  quit(save = "no", status = 0)
}

# ---- Dry run -----------------------------------------------------------------------
# The default. This is the only step in the pipeline that writes to a bucket other people
# read, so it plans by default and uploads only when asked, the same posture `make
# promote` takes toward production. The plan is written to the same log the real run
# writes, so the two are directly comparable.

if (!APPLY) {
  plan <- data.frame(
    gcs_path = to_upload$gcs_path, change_type = to_upload$change_type,
    old_version = to_upload$old_version, new_version = to_upload$new_version,
    would_remove = ifelse(to_upload$change_type %in% SUPERSEDING,
                          to_upload$staging_gcs_path, NA_character_),
    stringsAsFactors = FALSE
  )
  report_write(plan, file.path(LOG_DIR, paste0("upload_plan_", .timestamp(), ".tsv")))
  # Write the overwrite record here too, so the list of objects that would lose their old
  # bytes can be read BEFORE the run that loses them, not only afterwards.
  write_replace_readme(to_upload[to_upload$change_type == "REPLACE", , drop = FALSE],
                       file.path(LOG_DIR, "upload_replace_README.md"),
                       applied = FALSE, bucket = STAGING_BUCKET)
  message("\nDRY RUN — nothing was uploaded. Set APPLY=1 to write to ", STAGING_BUCKET, ".")
  quit(save = "no", status = 77)   # SKIP: the driver records it, the stage does not fail
}

# ---- Upload ------------------------------------------------------------------------

results <- vector("list", nrow(to_upload))

for (i in seq_len(nrow(to_upload))) {
  row <- to_upload[i, ]
  message("[", i, "/", nrow(to_upload), "] ", row$change_type, " ", row$rel_path)

  exit <- system(paste(GSUTIL, "cp", shQuote(row$local_path), shQuote(row$gcs_path)))

  # Post-upload verification. Read the object back and compare its md5 to the local file.
  stat_out <- system(paste(GSUTIL, "stat", shQuote(row$gcs_path)), intern = TRUE)
  md5_line <- grep("Hash (md5):", stat_out, fixed = TRUE, value = TRUE)
  remote_md5 <- if (length(md5_line)) md5_base64_to_hex(trimws(sub("^[^:]*:", "", md5_line[1])))
                else NA_character_
  verified <- !is.na(remote_md5) && identical(remote_md5, row$new_md5)

  if (!verified) {
    message("  MD5 MISMATCH")
    message("    local:  ", row$new_md5)
    message("    remote: ", remote_md5)
  }

  # Retire the superseded object. Only for a change type that moved the version, only after
  # the replacement verified, and never the object just written.
  removed <- NA
  if (row$change_type %in% SUPERSEDING && verified &&
      !is.na(row$staging_gcs_path) && !identical(row$staging_gcs_path, row$gcs_path)) {
    rm_exit <- system(paste(GSUTIL, "rm", shQuote(row$staging_gcs_path)))
    removed <- (rm_exit == 0)
    if (!removed) message("  WARNING: failed to remove superseded object: ", row$staging_gcs_path)
  }

  results[[i]] <- data.frame(
    gcs_path = row$gcs_path, superseded_gcs_path = row$staging_gcs_path,
    change_type = row$change_type, old_version = row$old_version,
    new_version = row$new_version, exit_code = exit,
    local_md5 = row$new_md5, remote_md5 = remote_md5,
    verified = verified, superseded_removed = removed,
    stringsAsFactors = FALSE
  )
}

all_results <- do.call(rbind, results)
log_path <- file.path(LOG_DIR, paste0("upload_log_", .timestamp(), ".tsv"))
report_write(all_results, log_path)

upload_failures  <- all_results[!all_results$verified, , drop = FALSE]
removal_failures <- all_results[all_results$change_type %in% SUPERSEDING &
                                !is.na(all_results$superseded_removed) &
                                !all_results$superseded_removed, , drop = FALSE]

report <- data.frame(
  status = c(if (nrow(upload_failures)) "FAIL" else "PASS",
             if (nrow(removal_failures)) "FAIL" else "PASS"),
  check  = c("upload:verified", "upload:superseded-removed"),
  detail = c(paste0(nrow(all_results) - nrow(upload_failures), "/", nrow(all_results),
                    " uploads md5-verified"),
             paste0(sum(all_results$superseded_removed, na.rm = TRUE), " superseded object(s) removed")),
  stringsAsFactors = FALSE
)
report_write(report, file.path(LOG_DIR, "upload_apply_report.tsv"))

# The overwrite record, written before the failure checks below so it survives a run that
# stops on a verification failure — a partially applied REPLACE set is exactly when knowing
# which objects were rewritten matters most.
replaced <- to_upload[to_upload$change_type == "REPLACE", , drop = FALSE]
write_replace_readme(replaced, file.path(LOG_DIR, "upload_replace_README.md"),
                     applied = TRUE, bucket = STAGING_BUCKET,
                     verified = all_results$verified[match(replaced$gcs_path,
                                                           all_results$gcs_path)])

if (nrow(upload_failures) > 0)
  stop(nrow(upload_failures), " file(s) failed MD5 verification. ",
       "The superseded versions were left in place. Inspect ", log_path)

if (nrow(removal_failures) > 0)
  stop(nrow(removal_failures), " superseded object(s) could not be removed — staging now ",
       "carries two versions of them. Inspect ", log_path)

# A verified apply makes diffs/latest.tsv stale by construction: every row it still lists
# as pending has just been uploaded, and nothing re-snapshots the bucket to say so. Step 05
# reads that same file to decide whether the bucket holds this release, so leaving it in
# place means step 05 sees the pre-upload pending count and refuses to validate — in a
# single APPLY=1 driver run it could never pass. Dropping it is what says "this diff no
# longer describes the bucket": step 05 validates when there is no diff to contradict it,
# and step 03 writes a fresh one on the next run. The timestamped diff_<ts>.tsv beside it
# is the durable record, so nothing is lost.
if (file.exists(diff_path)) {
  file.remove(diff_path)
  message("diffs/latest.tsv removed — it described the bucket before this upload.")
}

message("Step 04 complete: ", nrow(all_results), " uploads verified.")
