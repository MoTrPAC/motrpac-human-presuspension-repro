# Vendored helpers for the self-contained leaf-object builders.
#
# These are copied from MotrpacHumanPreSuspensionAnalysis so the preflight stage
# can build the dependency-free ("leaf", deps=-) data objects WITHOUT installing
# the Analysis package. Keep in sync with the upstream source:
#   .writeGMT                      -> internal package writer (GMT description = " ")
#   ome_available_list / .find_ome / .find_tissue  -> R/available_ome_tissue.R
#
# GMT line format (tab-separated), matching the committed *.gmt.gz files:
#   <set_name>\t \t<member1>\t<member2>...        (description field is a single space)

.writeGMT <- function(x, path) {
  stopifnot(!is.null(names(x)))
  lines <- vapply(seq_along(x), function(i) {
    paste(c(names(x)[i], " ", x[[i]]), collapse = "\t")
  }, character(1))
  # CRLF line endings to match the committed *.gmt.gz byte-for-byte (decompressed)
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeLines(lines, con, sep = "\r\n")
  invisible(path)
}

ome_available_list <- function() {
  # includes the two clinical omes prot-clinical / metab-t-clinical
  c(
    "prot-ol", "prot-ph", "prot-pr", "transcript-rna-seq",
    "epigen-methylcap-seq", "epigen-atac-seq",
    "metab-u-hilicpos", "metab-u-ionpneg", "metab-u-lrpneg", "metab-u-lrppos",
    "metab-u-rpneg", "metab-u-rppos", "metab-t-amines", "metab-t-conv",
    "metab-t-imm-crt", "metab-t-oxylipneg", "metab-t-tca", "metab-t-nuc",
    "metab-t-acoa", "metab-t-ka",
    "prot-clinical", "metab-t-clinical"
  )
}

.find_ome <- function(file_path) {
  for (ome in ome_available_list())
    if (grepl(ome, file_path)) return(ome)
  NULL
}

.find_tissue <- function(file_path) {
  tissue_combinations <- list(
    "t02" = "blood", "t03" = "blood", "t04" = "blood", "t05" = "blood",
    "t06" = "muscle", "t10" = "muscle", "t07" = "adipose", "t11" = "adipose"
  )
  for (code in names(tissue_combinations))
    if (grepl(code, file_path)) return(tissue_combinations[[code]])
  NULL
}

# --- shared paths for the builders ------------------------------------------
# build_data_objects.sh exports PREFLIGHT_DR (absolute data-raw dir),
# PREFLIGHT_DATA_DIR (rds-only, R-package-standard) and PREFLIGHT_OUTPUT_DIR
# (everything else: gmt.gz, json). Builders resolve through these so they run
# regardless of the caller's working directory.
#   .out_dir()    -> data/    (.rds data objects ONLY)
#   .output_dir() -> output/  (gmt.gz, json, other non-rds artifacts)
.dataraw_dir <- function() {
  d <- Sys.getenv("PREFLIGHT_DR", "")
  if (!nzchar(d)) d <- normalizePath(".")
  d
}
.out_dir <- function() {
  d <- Sys.getenv("PREFLIGHT_DATA_DIR", "")
  if (!nzchar(d)) d <- file.path(dirname(.dataraw_dir()), "data")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  normalizePath(d)
}
.output_dir <- function() {
  d <- Sys.getenv("PREFLIGHT_OUTPUT_DIR", "")
  if (!nzchar(d)) d <- file.path(dirname(.dataraw_dir()), "output")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  normalizePath(d)
}
# Compress a written .gmt to .gmt.gz in the OUTPUT dir (overwrite-safe).
.gmt_out <- function(named_list, gmt_basename) {
  out <- .output_dir()
  gmt <- file.path(out, gmt_basename)
  .writeGMT(named_list, gmt)
  gz <- paste0(gmt, ".gz")
  if (file.exists(gz)) file.remove(gz)
  R.utils::gzip(gmt, destname = gz, remove = TRUE)
  gz
}
