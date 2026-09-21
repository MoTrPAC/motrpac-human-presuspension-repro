#!/usr/bin/env Rscript
# Leaf data object: QUANTID_BUCKET_FILES  (gated — reads the quant-id bucket)
#
# Catalogs the raw post-quantification inputs the qc-norm stems consume. Lists the
# quant-id bucket (QUANTID_BUCKET in config/pipeline.env) once and parses every data
# file's gs:// path into the tags that quantid_path() (lib/qc_helpers.R) filters on:
#   assay_class, tissue_code, data_details, content_type, gcs_path, is_qa_qc
# Writes the catalog as BOTH quantid_bucket_files.json (the human-readable source) and
# QUANTID_BUCKET_FILES.rds (loaded at run time). Regenerable from the bucket — no
# longer a hand-vendored file.
#
# Path layout:  <QUANTID_BUCKET>/<version>/<group>/<tissue_code|qa-qc>/<ome>/<file>
# Filename:     motrpac_human-precovid_[<tissue_code>_]<ome>_<data_details>_v<ver>.<ext>
suppressWarnings(suppressMessages({ library(jsonlite) }))

out_dir  <- Sys.getenv("PREFLIGHT_OUTPUT_DIR")
data_dir <- Sys.getenv("PREFLIGHT_DATA_DIR")
gsutil   <- Sys.getenv("GSUTIL", "gsutil")
stopifnot(nzchar(out_dir), nzchar(data_dir))

# QUANTID_BUCKET from config (Sys.getenv wins if a driver exported it).
quantid_bucket <- Sys.getenv("QUANTID_BUCKET")
if (!nzchar(quantid_bucket)) {
  root <- dirname(dirname(dirname(out_dir)))          # .../motrpac-human-presuspension-repro
  env  <- file.path(root, "config", "pipeline.env")
  ln   <- grep("^(export )?QUANTID_BUCKET=", readLines(env, warn = FALSE), value = TRUE)
  if (length(ln)) quantid_bucket <- gsub('^"|"$', "", sub(".*QUANTID_BUCKET=", "", ln[length(ln)]))
}
if (!nzchar(quantid_bucket)) stop("QUANTID_BUCKET not set and not found in config/pipeline.env")

# ---- list the bucket + keep data files -------------------------------------
all <- system(paste(gsutil, "ls -R", shQuote(quantid_bucket)), intern = TRUE)
all <- all[grepl("\\.(txt|csv)(\\.gz)?$", all)]
if (!length(all)) stop("no data files found under ", quantid_bucket)

# ---- parse each path into the catalog tags ---------------------------------
parse_one <- function(gcs) {
  fn <- basename(gcs)
  after_ver <- sub(".*/v[0-9.]+/", "", gcs)                          # segments below v<ver>/
  dirs <- head(strsplit(after_ver, "/", fixed = TRUE)[[1]], -1L)     # drop the filename
  is_qa_qc <- "qa-qc" %in% dirs
  tissue_code <- { tc <- dirs[grepl("^t[0-9]+-", dirs)]; if (length(tc)) tc[1] else NA_character_ }
  core <- sub("^motrpac_human-precovid_", "", fn)
  core <- sub("_v[0-9.]+(-[a-z0-9]+)?\\.(txt|csv)(\\.gz)?$", "", core)   # strip version + ext
  if (is_qa_qc) {                                                    # <ome>_qa-qc-metrics
    assay_class  <- sub("_qa-qc-metrics$", "", core)
    data_details <- "qa-qc-metrics"
  } else {
    assay_class  <- tail(dirs, 1L)                                   # the ome directory
    prefix <- paste0(if (!is.na(tissue_code)) paste0(tissue_code, "_") else "", assay_class, "_")
    data_details <- sub(paste0("^", prefix), "", core)
  }
  data.frame(assay_class = assay_class, tissue_code = tissue_code, data_details = data_details,
             content_type = NA_character_, gcs_path = gcs, is_qa_qc = is_qa_qc,
             stringsAsFactors = FALSE)
}
files <- do.call(rbind, lapply(all, parse_one))
files <- files[order(files$assay_class, files$tissue_code, files$data_details), ]
rownames(files) <- NULL

# sanity: the (assay_class, tissue_code, data_details) key must be unique so
# quantid_path() (which does not filter by content_type) resolves to one file.
k <- files[c("assay_class", "tissue_code", "data_details")]
dup <- files[duplicated(k) | duplicated(k, fromLast = TRUE), ]
if (nrow(dup))
  warning(sprintf("%d non-unique (assay_class,tissue_code,data_details) rows — quantid_path may be ambiguous:\n%s",
                  nrow(dup), paste(utils::head(dup$gcs_path, 6), collapse = "\n")))

# ---- write json (source of truth) + rds (run-time object) ------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
jsonlite::write_json(list(source_bucket = quantid_bucket, files = files),
                     file.path(out_dir, "quantid_bucket_files.json"),
                     pretty = TRUE, auto_unbox = TRUE, na = "null")
saveRDS(files, file.path(data_dir, "QUANTID_BUCKET_FILES.rds"))
message("QUANTID_BUCKET_FILES: ", nrow(files), " files cataloged from ", quantid_bucket,
        " -> quantid_bucket_files.json + QUANTID_BUCKET_FILES.rds")
