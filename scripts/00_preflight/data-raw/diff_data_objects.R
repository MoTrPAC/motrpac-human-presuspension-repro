#!/usr/bin/env Rscript
# Diff each built data object against the version in the package that creates it.
# .rds objects  -> <pkg>/data/<NAME>.rda        (resolved by where the object lives)
# .gmt.gz files -> <analysis>/data-raw/gmt_processing/gmt_files/<name>.gmt.gz (membership)
# Objects with no package version (e.g. the repo-local quant-id catalog) are noted.
# Env: PREFLIGHT_DATA_DIR, PREFLIGHT_OUTPUT_DIR, ANALYSIS_PKG_REPO, DATA_PKG_REPO.

data_dir   <- Sys.getenv("PREFLIGHT_DATA_DIR")
output_dir <- Sys.getenv("PREFLIGHT_OUTPUT_DIR")
ANALYSIS   <- Sys.getenv("ANALYSIS_PKG_REPO")
pkgs       <- Filter(nzchar, c(ANALYSIS, Sys.getenv("DATA_PKG_REPO")))

pkg_label <- function(path) if (grepl("Analysis", path)) "analysis" else
                            if (grepl("Data", path)) "data" else basename(path)

find_rda <- function(name) {
  for (p in pkgs) { f <- file.path(p, "data", paste0(name, ".rda"))
    if (file.exists(f)) return(f) }
  NA_character_
}
load_rda <- function(rda, name) {
  e <- new.env(); load(rda, envir = e)
  if (name %in% ls(e)) e[[name]] else e[[ls(e)[1]]]
}

cmp <- function(a, b) {
  if (identical(a, b)) return("identical")
  if (is.data.frame(a) && is.data.frame(b)) {
    if (!identical(dim(a), dim(b)))
      return(sprintf("DIFFER (dim %s vs %s)",
                     paste(dim(a), collapse = "x"), paste(dim(b), collapse = "x")))
    ord <- function(d) { d <- d[do.call(order, lapply(d, as.character)), , drop = FALSE]
                         rownames(d) <- NULL; d }
    if (isTRUE(all.equal(ord(a), ord(b), check.attributes = FALSE)))
      return("identical (row order differs)")
    return("DIFFER (content)")
  }
  r <- all.equal(a, b)
  if (isTRUE(r)) return("near-equal (all.equal)") else return(sprintf("DIFFER (%d)", length(r)))
}

rows <- list(); add <- function(...) rows[[length(rows) + 1]] <<- c(...)

# ---- .rds data objects -----------------------------------------------------
for (f in sort(list.files(data_dir, pattern = "\\.rds$", full.names = TRUE))) {
  name <- sub("\\.rds$", "", basename(f))
  rda  <- find_rda(name)
  if (is.na(rda)) { add(name, "-", "no package version"); next }
  res <- tryCatch(cmp(readRDS(f), load_rda(rda, name)),
                  error = function(e) paste("ERROR:", conditionMessage(e)))
  add(name, pkg_label(rda), res)
}

# ---- .gmt.gz (order-independent set membership) ----------------------------
gdir <- file.path(ANALYSIS, "data-raw", "gmt_processing", "gmt_files")
have_tmsig <- requireNamespace("TMSig", quietly = TRUE)
for (f in sort(list.files(output_dir, pattern = "\\.gmt\\.gz$", full.names = TRUE))) {
  base <- basename(f); pk <- file.path(gdir, base)
  if (!file.exists(pk)) { add(base, "-", "no package version"); next }
  if (!have_tmsig)      { add(base, "analysis", "skip (TMSig missing)"); next }
  a <- TMSig::readGMT(f, check = FALSE); b <- TMSig::readGMT(pk, check = FALSE)
  if (!setequal(names(a), names(b))) { add(base, "analysis", "DIFFER (set names)"); next }
  mem_ok <- all(vapply(names(a), function(n) setequal(a[[n]], b[[n]]), logical(1)))
  add(base, "analysis", if (mem_ok) "identical (membership)" else "DIFFER (membership)")
}

# ---- non-.rds/.gmt outputs with no package version -------------------------
jf <- file.path(output_dir, "quantid_bucket_files.json")
if (file.exists(jf)) add("quantid_bucket_files.json", "-", "no package version")

cat(sprintf("%-32s %-9s %s\n", "OBJECT", "PACKAGE", "DIFF"))
cat(sprintf("%-32s %-9s %s\n", "------", "-------", "----"))
for (r in rows) cat(sprintf("%-32s %-9s %s\n", r[1], r[2], r[3]))
