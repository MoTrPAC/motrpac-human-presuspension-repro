#!/usr/bin/env Rscript
# Step 01 — timestamped MD5 inventory of the production and staging buckets.
#
# Vendored from google_cloud_bucket_checks/01_snapshot_production.R. Two adaptations,
# both explained in the step README: the staging bucket is inventoried alongside
# production, and the listing is a single `gsutil ls -L` pass rather than one
# `gsutil stat` per object.

source(file.path(Sys.getenv("STAGE_LIB"), "bucket_helpers.R"))

ts <- .timestamp()

snapshot_one <- function(label, bucket) {
  message("Listing ", label, " bucket: ", bucket)
  inv <- gcs_inventory(bucket)
  if (nrow(inv) == 0)
    stop("No files found in ", label, " bucket (", bucket, "). ",
         "Check the path in config/pipeline.env and your gcloud credentials.")

  inv$rel_path           <- sub(paste0("^", bucket, "/"), "", inv$gcs_path)
  inv$snapshot_timestamp <- ts
  inv$bucket             <- bucket

  missing_md5 <- sum(is.na(inv$md5))
  if (missing_md5 > 0)
    warning(missing_md5, " object(s) in ", label, " have no md5 hash — ",
            "composite objects cannot be content-compared and will read as MODIFIED")

  path <- file.path(SNAPSHOT_DIR, paste0("snapshot_", label, "_", ts, ".tsv"))
  report_write(inv, path)
  file.copy(path, file.path(SNAPSHOT_DIR, paste0("latest_", label, ".tsv")), overwrite = TRUE)

  message("  ", nrow(inv), " files inventoried")
  inv
}

prod    <- snapshot_one("production", PRODUCTION_BUCKET)
staging <- snapshot_one("staging",    STAGING_BUCKET)

report <- data.frame(
  status = "PASS",
  check  = c("snapshot:production", "snapshot:staging"),
  detail = c(paste0(nrow(prod), " files at ", PRODUCTION_BUCKET),
             paste0(nrow(staging), " files at ", STAGING_BUCKET)),
  stringsAsFactors = FALSE
)
report_write(report, file.path(LOG_DIR, "upload_snapshot_report.tsv"))

message("Step 01 complete: production ", nrow(prod), ", staging ", nrow(staging), " files.")
