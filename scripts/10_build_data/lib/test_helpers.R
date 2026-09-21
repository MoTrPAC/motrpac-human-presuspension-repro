# Lightweight test harness + assertion catalog for the 10_build_data test suite.
# No testthat. A step's <name>_tests.R sources this, calls assert()/check()/skip(),
# then finish() to print a summary and exit non-zero on any FAIL.
#
# All three tiers are available; a check whose prerequisite is genuinely absent
# (package .rda not found, no gsutil / bucket access) records SKIP, not FAIL.
suppressWarnings(suppressMessages({ library(data.table) }))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ── paths / shared locations ──────────────────────────────────────────────────
.T_ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.T_ROOT)) {
  .T_ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(.T_ROOT, "config", "pipeline.env")) && dirname(.T_ROOT) != .T_ROOT)
    .T_ROOT <- dirname(.T_ROOT)
}
# From a linked worktree dirname(.T_ROOT) is .claude/worktrees/, not the GitHub dir, so the
# package checkouts below would not be found. The main checkout is the parent of the shared
# git dir; same resolution config/pipeline.env uses for GITHUB_ROOT.
.T_GH        <- dirname(.T_ROOT)
if (!dir.exists(file.path(.T_GH, "MotrpacHumanPreSuspensionAnalysis"))) {
  .common <- suppressWarnings(tryCatch(
    system2("git", c("-C", shQuote(.T_ROOT), "rev-parse", "--path-format=absolute",
                     "--git-common-dir"), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)))
  if (length(.common) == 1L && nzchar(.common) && dir.exists(.common))
    .T_GH <- dirname(dirname(.common))
}
.T_DATA      <- file.path(.T_ROOT, "scripts", "10_build_data", "data")
.T_FREEZE    <- file.path(.T_ROOT, "staging", "freeze")
.T_RAW       <- file.path(.T_ROOT, "staging", "raw-files")   # persistent bucket-download cache
.T_PREFLIGHT <- file.path(.T_ROOT, "scripts", "00_preflight", "data")
.T_ANALYSIS  <- file.path(.T_GH, "MotrpacHumanPreSuspensionAnalysis")
.T_DATAPKG   <- file.path(.T_GH, "MotrpacHumanPreSuspensionData")
.T_PKGS      <- Filter(dir.exists, c(.T_ANALYSIS, .T_DATAPKG))

# Read a value from config/pipeline.env (Sys.getenv wins if a driver exported it).
.pipeline_env_val <- function(name, default = "") {
  v <- Sys.getenv(name); if (nzchar(v)) return(v)
  f <- file.path(.T_ROOT, "config", "pipeline.env"); if (!file.exists(f)) return(default)
  ln <- grep(paste0("(^|export )", name, "="), readLines(f, warn = FALSE), value = TRUE)
  if (!length(ln)) return(default)
  gsub('^"|"$', "", sub(paste0(".*", name, "="), "", ln[length(ln)]))
}
# Reference freeze bucket for the diff-vs-bucket tier — from STAGING_BUCKET (config).
.T_BUCKET <- .pipeline_env_val("STAGING_BUCKET")

# ── harness ───────────────────────────────────────────────────────────────────
.TH <- new.env(); .TH$rows <- list()
.record <- function(status, desc, detail = "") {
  .TH$rows[[length(.TH$rows) + 1]] <- list(status = status, desc = desc, detail = detail)
  tag <- switch(status, PASS = "[ OK ]", FAIL = "[FAIL]", SKIP = "[SKIP]",
                INFO = "[INFO]", WARN = "[WARN]")
  cat(sprintf("  %s %s%s\n", tag, desc, if (nzchar(detail)) paste0(" — ", detail) else ""))
  invisible(status)
}
assert <- function(desc, cond, fail_detail = "") {
  .record(if (isTRUE(cond)) "PASS" else "FAIL", desc, if (isTRUE(cond)) "" else fail_detail)
}
check <- function(desc, expr) {
  res <- tryCatch(isTRUE(expr), error = function(e) structure(FALSE, msg = conditionMessage(e)))
  if (isTRUE(res)) .record("PASS", desc) else .record("FAIL", desc, attr(res, "msg") %||% "condition FALSE")
}
skip   <- function(desc, reason = "") .record("SKIP", desc, reason)
# report(): log an informational finding (e.g. a diff-vs-package result) that must
# NOT fail the step — it never counts toward the exit status.
report <- function(desc, detail = "") .record("INFO", desc, detail)
# warn(): non-fatal, louder than report(). For findings the operator must see but which must
# not stop the build — a truncated DA table is real, but step 09 is mid-migration and several
# stems have not been rerun, so failing here would block every build for a known state.
warn <- function(desc, detail = "") .record("WARN", desc, detail)
finish <- function() {
  st <- vapply(.TH$rows, `[[`, "", "status")
  cat(sprintf("SUMMARY: %d PASS, %d FAIL, %d WARN, %d SKIP, %d INFO\n",
              sum(st == "PASS"), sum(st == "FAIL"), sum(st == "WARN"),
              sum(st == "SKIP"), sum(st == "INFO")))
  quit(save = "no", status = if (sum(st == "FAIL") > 0) 1L else 0L)
}

# ── local object loaders ──────────────────────────────────────────────────────
load_built <- function(name) {
  f <- file.path(.T_DATA, paste0(name, ".rda"))
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e); e[[name]] %||% e[[ls(e)[1]]]
}
load_preflight <- function(name) {
  f <- file.path(.T_PREFLIGHT, paste0(name, ".rds")); if (!file.exists(f)) NULL else readRDS(f)
}
# For the nested [[tissue]][[ome]] list, source lib/qc_helpers.R and call load_qc_local()
# directly — see qc_object_tests.R. No wrapper here: it is the same function the DA stems
# call, so the tests should reach it the same way, with its own defaults.

# ── structural assertions ─────────────────────────────────────────────────────
assert_cols <- function(obj_desc, df, required, classes = NULL) {
  miss <- setdiff(required, colnames(df))
  assert(sprintf("%s: columns present", obj_desc), length(miss) == 0,
         sprintf("missing: %s", paste(miss, collapse = ", ")))
  for (cc in names(classes)) if (cc %in% colnames(df)) {
    assert(sprintf("%s: %s is %s", obj_desc, cc, classes[[cc]]),
           classes[[cc]] %in% class(df[[cc]]), sprintf("got %s", class(df[[cc]])[1]))
  }
}
assert_levels <- function(desc, x, expected) {
  assert(desc, identical(levels(x), expected), sprintf("got [%s]", paste(levels(x), collapse = ", ")))
}
assert_subset <- function(desc, a, b) {  # a ⊆ b
  extra <- setdiff(a, b)
  assert(desc, length(extra) == 0, sprintf("%d not in target, e.g. %s", length(extra),
         paste(utils::head(extra, 5), collapse = ", ")))
}
assert_unique <- function(desc, x) {
  d <- x[duplicated(x)]
  assert(desc, length(d) == 0, sprintf("%d dup(s), e.g. %s", length(d), paste(utils::head(unique(d), 5), collapse = ", ")))
}
assert_na <- function(desc, x, max_frac = 1) {
  fr <- mean(is.na(x)); assert(desc, fr <= max_frac, sprintf("NA frac %.3f > %.3f", fr, max_frac))
}
assert_completeness <- function(desc, present, expected) {
  miss <- setdiff(expected, present)
  assert(desc, length(miss) == 0, sprintf("missing %d: %s", length(miss),
         paste(utils::head(miss, 8), collapse = ", ")))
}
# no two freeze files share a (token,type) differing only by _vX.Y version
assert_no_dup_versions <- function(freeze_dir = .T_FREEZE) {
  files <- list.files(freeze_dir, pattern = "\\.txt$", recursive = TRUE)
  key <- sub("_v[0-9.]+\\.txt$", "", basename(files))
  dups <- unique(key[duplicated(key)])
  assert("no duplicate-version freeze files", length(dups) == 0,
         sprintf("%d dup key(s), e.g. %s", length(dups), paste(utils::head(dups, 5), collapse = ", ")))
}

# the ome x tissue combos that should exist (qc / DA), from OME_TISSUE_CODE
expected_tissue_omes <- function(include_clinical = TRUE) {
  otc <- load_preflight("OME_TISSUE_CODE")
  if (is.null(otc)) return(NULL)
  keep <- !grepl("^epigen-|^lab-", otc$ome)
  if (!include_clinical) keep <- keep & !grepl("clinical", otc$ome)
  otc <- otc[keep, ]
  paste0(otc$tissue, ".", otc$ome)
}

# per-ome required covariate columns (drives process_covariates)
covariate_cols_for <- function(ome, tissue) {
  csv <- file.path(.T_ROOT, "scripts", "00_preflight", "data-raw", "sources", "covariates_pre_cawg.csv")
  if (!file.exists(csv)) return(character(0))
  cv <- utils::read.csv(csv, stringsAsFactors = FALSE)
  cv <- cv[cv$ome == ome & (cv$tissue == "all" | cv$tissue == tissue), , drop = FALSE]
  unique(as.character(cv$covariate))
}

# validate one *_QC leaf against the load_qc / DA contract
QC_SAMPLE_REQUIRED <- c("vialLabel", "randomGroupCode", "Sex", "Timepoint",
                        "Age_3_groups", "visitcode", "pid")
TIMEPOINT_LEVELS <- c("pre_exercise", "during_20_min", "during_40_min", "post_10_min",
                      "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
assert_qc_leaf <- function(tissue, ome, leaf) {
  d <- sprintf("%s/%s", tissue, ome)
  assert(sprintf("%s: has qc_norm+sample+feature_metadata", d),
         all(c("qc_norm", "sample_metadata", "feature_metadata") %in% names(leaf)),
         sprintf("has [%s]", paste(names(leaf), collapse = ", ")))
  qn <- leaf$qc_norm; sm <- leaf$sample_metadata
  if (is.null(qn) || is.null(sm)) return(invisible())
  assert(sprintf("%s: qc_norm is data.frame", d), is.data.frame(qn))
  # A sample column that is entirely NA carries no type and so is logical, not
  # numeric; the published package objects contain such columns, so only real
  # character/factor contamination counts as a defect here.
  nonnum <- names(qn)[!vapply(qn, function(x) is.numeric(x) || all(is.na(x)), logical(1))]
  assert(sprintf("%s: qc_norm values numeric", d), length(nonnum) == 0,
         sprintf("non-numeric column(s): %s", paste(utils::head(nonnum, 5), collapse = ", ")))
  assert_unique(sprintf("%s: qc_norm feature rownames unique", d), rownames(qn))
  # sample_metadata required + per-ome covariate columns
  cov <- covariate_cols_for(ome, tissue)
  assert_cols(d, sm, unique(c(QC_SAMPLE_REQUIRED, cov)))
  assert(sprintf("%s: visitcode has ADU_BAS", d), "ADU_BAS" %in% as.character(sm$visitcode))
  if (is.factor(sm$Timepoint))
    assert_subset(sprintf("%s: Timepoint levels valid", d), levels(sm$Timepoint), TIMEPOINT_LEVELS)
  # alignment: qc_norm colnames ⊆ vialLabel ; feature_metadata id ↔ rownames
  assert_subset(sprintf("%s: qc_norm cols ⊆ vialLabel", d),
                colnames(qn), as.character(sm$vialLabel))
  fm <- leaf$feature_metadata
  if (!is.null(fm) && "id" %in% colnames(fm))
    assert(sprintf("%s: feature_metadata$id ⊇ qc_norm rows", d),
           all(rownames(qn) %in% as.character(fm$id)), "some rownames not in feature_metadata$id")
}

# ── diff vs package .rda (vendored from 00_preflight/data-raw/diff_data_objects.R) ─
.find_rda <- function(name) {
  for (p in .T_PKGS) { f <- file.path(p, "data", paste0(name, ".rda")); if (file.exists(f)) return(f) }
  NA_character_
}
.load_rda <- function(rda, name) { e <- new.env(); load(rda, envir = e); e[[name]] %||% e[[ls(e)[1]]] }
.cmp <- function(a, b) {
  if (identical(a, b)) return("identical")
  if (is.data.frame(a) && is.data.frame(b)) {
    if (!identical(dim(a), dim(b)))
      return(sprintf("DIFFER (dim %s vs %s)", paste(dim(a), collapse = "x"), paste(dim(b), collapse = "x")))
    ord <- function(d) { d <- d[do.call(order, lapply(d, as.character)), , drop = FALSE]; rownames(d) <- NULL; d }
    if (isTRUE(all.equal(ord(a), ord(b), check.attributes = FALSE))) return("identical (row order differs)")
    return("DIFFER (content)")
  }
  r <- all.equal(a, b); if (isTRUE(r)) "near-equal (all.equal)" else sprintf("DIFFER (%d)", length(r))
}
# diff vs the shipped package object. This is a REPORT (INFO), never a FAIL — a
# reproduction difference is logged for review, not a build-stopping error.
# (`tolerant` kept for call-site compatibility; ignored.)
diff_vs_package <- function(name, obj = NULL, tolerant = FALSE) {
  rda <- .find_rda(name)
  if (is.na(rda)) return(report(sprintf("diff vs package: %s", name), "no package .rda"))
  obj <- obj %||% load_built(name)
  if (is.null(obj)) return(report(sprintf("diff vs package: %s", name), "not built locally"))
  res <- tryCatch(.cmp(obj, .load_rda(rda, name)), error = function(e) paste("ERROR:", conditionMessage(e)))
  report(sprintf("diff vs package: %s", name), res)
}

# ── diff vs GCS bucket freeze (vendored per_feat_cor + gsutil align) ───────────
.have_gsutil <- function() nzchar(Sys.which("gsutil"))
.per_feat_cor <- function(A, B) {
  n <- nrow(A)
  if (n > 40000 && !anyNA(A) && !anyNA(B)) {
    A0 <- A - rowMeans(A); B0 <- B - rowMeans(B)
    r <- rowSums(A0 * B0) / sqrt(rowSums(A0^2) * rowSums(B0^2)); r[!is.finite(r)] <- NA_real_; return(r)
  }
  vapply(seq_len(n), function(i) {
    x <- A[i, ]; y <- B[i, ]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) NA_real_ else cor(x[ok], y[ok])
  }, numeric(1))
}
# Cached one-shot listing of the reference bucket (avoids a per-file gsutil ls).
.BUCKET_CACHE <- new.env()
.bucket_files <- function() {
  if (is.null(.BUCKET_CACHE$f)) {
    if (!.have_gsutil() || !nzchar(.T_BUCKET)) { .BUCKET_CACHE$f <- character(0) } else {
      x <- tryCatch(system(sprintf("gsutil ls -R '%s'", .T_BUCKET), intern = TRUE, ignore.stderr = TRUE),
                    error = function(e) character(0))
      .BUCKET_CACHE$f <- x[grepl("\\.txt$", x)]
    }
  }
  .BUCKET_CACHE$f
}
# Compare one local qc-norm freeze matrix to its counterpart on STAGING_BUCKET by
# per-feature correlation across shared features x samples. Reports the median cor in
# the detail; FAILs only below min_cor (a GROSS-misalignment gate — the known
# re-referencing/version differences, e.g. transcript ~0.92, are acceptable, so this
# catches catastrophic errors like a mis-aligned cbind, not benign drift).
diff_freeze_vs_bucket <- function(local, min_cor = 0.5) {
  desc <- sprintf("diff vs bucket: %s", basename(local))
  remote_all <- .bucket_files()
  if (length(remote_all) == 0) return(skip(desc, "no bucket listing (gsutil/bucket unavailable)"))
  rel <- sub(paste0("^", .T_FREEZE, "/"), "", local)
  keyb <- sub("_v[0-9.]+\\.txt$", "", basename(local))
  cand <- remote_all[sub("_v[0-9.]+\\.txt$", "", basename(remote_all)) == keyb &
                     grepl(paste0("/", dirname(rel), "/"), remote_all)]
  if (length(cand) != 1) return(skip(desc, sprintf("%d bucket matches on %s", length(cand), .T_BUCKET)))
  # download the reference file ONCE into staging/raw-files and reuse it on later runs
  dir.create(.T_RAW, showWarnings = FALSE, recursive = TRUE)
  cached <- file.path(.T_RAW, basename(cand[1]))
  if (!file.exists(cached) &&
      system(sprintf("gsutil cp '%s' '%s'", cand[1], cached), ignore.stdout = TRUE, ignore.stderr = TRUE) != 0)
    return(skip(desc, "gsutil cp failed"))
  # fread, not read.csv: this diff loads two full qc-norm matrices at once.
  o <- fread(local,  sep = "\t", header = TRUE, data.table = FALSE)
  s <- fread(cached, sep = "\t", header = TRUE, data.table = FALSE)
  cf <- intersect(o[[1]], s[[1]]); cs <- intersect(colnames(o)[-1], colnames(s)[-1])
  if (length(cf) < 3 || length(cs) < 3) return(skip(desc, "too few shared features/samples"))
  med <- stats::median(.per_feat_cor(as.matrix(o[match(cf, o[[1]]), cs, drop = FALSE]),
                                     as.matrix(s[match(cf, s[[1]]), cs, drop = FALSE])), na.rm = TRUE)
  assert(sprintf("%s (median per-feat cor %.3f)", desc, med), is.finite(med) && med >= min_cor,
         sprintf("median cor %.3f < %.2f — likely misalignment", med, min_cor))
}
