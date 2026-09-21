#!/usr/bin/env Rscript
# Tests for the DA FREEZE (step 09). Downstream: step 11 (*_SUM_STATS) reads these txt files
# to decide which features may be summarised, and the staging-bucket upload ships them as-is.
#
# This step ran without tests for the whole first build cycle, and it cost correctness. Under
# a broken variancePartition/lme4 stack, dream() silently dropped every feature whose
# expression carried NAs: muscle prot-ph shipped DA for 1,122 of 17,955 features and muscle
# prot-pr 3,404 of 6,197, with nothing logged anywhere — the row counts simply came out
# lower. A second defect, the BiocParallel SOCK backend losing one feature per worker, hid
# the same way. Section 5 exists to make that class of silent truncation visible.
#
# SEVERITY SPLIT. Sections 2-4 are hard assertions: they describe a malformed file, which is
# always a defect. Sections 5 and 6 only warn or report, and never fail the build — they
# compare against a moving reference while step 09 is mid-migration, and several stems are
# knowingly red until they are rerun.
#
# Ported from data-raw/google_cloud_bucket_checks/, whose DA schemas have never actually run
# against a DA file: its manifest sets da_details = NA for every metab row and drops them,
# and its README says DA support "will come later". The corrections that port needed are
# recorded at each section.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
suppressWarnings(suppressMessages({ library(data.table); library(dplyr) }))
cat("== DA freeze tests ==\n")

# fread, not read.csv, in this file: sections 4-7 need a header, one column, or a subset of
# columns from tables that reach 1.2 GB, and read.csv cannot read part of a file. The same
# documented exception qc_norm_freeze_tests.R takes, for the same reason.

# Row-level checks are skipped above this size; a bare header/one-column read still runs.
MAX_MB <- suppressWarnings(as.numeric(Sys.getenv("DA_TESTS_MAX_MB", "512")))
if (!is.finite(MAX_MB)) MAX_MB <- 512

# ── 1. discovery ──────────────────────────────────────────────────────────────
# Two data_details strings, not one: the methylcap DA is a vendored malax-glmm table, so a
# glob on "_da_dream-" alone silently drops three of the largest files in the freeze.
DA_SUFFIX <- "_da_(dream|malax-glmm)-acute_v[0-9.]+\\.txt$"

assert_no_dup_versions()

da_files <- list.files(.T_FREEZE, pattern = DA_SUFFIX, recursive = TRUE, full.names = TRUE)
assert("DA freeze files present", length(da_files) > 0, "staging/freeze/*/da/ empty — run step 09 first")
if (length(da_files) == 0) finish()

token_of <- function(p) sub(DA_SUFFIX, "", sub("^human-precovid-sed-adu_", "", basename(p)))
ome_of   <- function(p) sub("^[^_]+_", "", token_of(p))
version_of <- function(p) numeric_version(sub("^.*_v([0-9.]+)\\.txt$", "\\1", basename(p)))

# One file per token for the token-keyed sections. assert_no_dup_versions() above already
# fails on a stale duplicate, so this is a fallback that keeps the rest of the suite running
# on the newer file rather than an arbitrary one.
all_tokens <- token_of(da_files)
newest <- vapply(split(seq_along(da_files), all_tokens),
                 function(ix) ix[which.max(xtfrm(version_of(da_files[ix])))], integer(1))
da_by_token <- stats::setNames(da_files[newest], names(newest))
report("DA freeze inventory",
       sprintf("%d file(s), %d token(s), %d over %d MB",
               length(da_files), length(da_by_token),
               sum(file.size(da_files) / 1048576 > MAX_MB), MAX_MB))

# ── 2. column schemas ─────────────────────────────────────────────────────────
# Source of truth is .convert_dream_output()'s closing select(), not upstream's vectors.
# Corrections that port needed, all of them upstream bugs rather than pipeline divergence:
#   - assay, z.std and AveExpr added. Upstream's DA_COLS_CORE omits all three although every
#     stem emits them, so its schema would pass a file that had lost the assay label.
#   - the separate prot / transcript / clinical schemas collapse to one vector, because one
#     select() produces all of them and they cannot differ.
#   - platform KEPT for metab. Upstream requires it and this pipeline never emitted it, which
#     is why generate_metab_da.R now does rather than the schema being weakened to match.
#   - epigen-atac-seq covered at all; upstream has no key for it.
# One interval pair, directly after logFC. topTable's own CI.L/CI.R stay absent, because
# variancePartition derives every contrast's bound from the first contrast's df;
# run_dream() passes confint = FALSE and rebuilds the interval itself, subsetting df.total
# by feature AND contrast. The corrected pair keeps the _calculated suffix so the bare
# names remain forbidden downstream. See .convert_dream_output() in da_common.R.
DA_COLS_DREAM  <- c("assay", "feature_id", "z.std", "logFC",
                    "CI.L_calculated", "CI.R_calculated",
                    "degrees_of_freedom", "logLik", "t", "AveExpr",
                    "p_value", "adj_p_value", "contrast", "full_model")
# generate_metab_da.R and generate_clinical_da.R relocate platform directly after assay, so
# this is a position, not an append — the order below is the contract, not just the column set.
DA_COLS_METAB  <- append(DA_COLS_DREAM, "platform", after = 1L)
DA_COLS_METHYL <- c("assay", "feature_id", "methylation_diff", "t", "AveExpr",
                    "p_value", "adj_p_value", "contrast", "full_model")

# Every metab DA file carries assay = "metab" with the ome in a platform column, including
# metab-t-clinical, which generate_clinical_da.R writes rather than the metab stem.
is_metab_platform <- function(ome) grepl("^metab", ome)
# The subset the metab stem fits. metab-t-clinical is not one of them: it is fit from the
# clinical *_QC objects, which are loaded whole rather than through the metab filters.
is_metab_stem_ome <- function(ome) is_metab_platform(ome) && ome != "metab-t-clinical"
schema_for <- function(ome) {
  if (ome == "epigen-methylcap-seq") return(DA_COLS_METHYL)
  if (is_metab_platform(ome))        return(DA_COLS_METAB)
  return(DA_COLS_DREAM)
}
header_of <- function(f) names(fread(f, nrows = 0, sep = "\t", header = TRUE))

for (f in da_files) {
  b <- basename(f); ome <- ome_of(f); want <- schema_for(ome); got <- header_of(f)
  if (ome == "epigen-methylcap-seq") {
    # vendored artifact rather than something this pipeline writes, so extra columns are
    # reported rather than treated as a schema break.
    assert(sprintf("DA %s: required columns", b), all(want %in% got),
           sprintf("missing: %s", paste(setdiff(want, got), collapse = ", ")))
    if (length(setdiff(got, want)))
      report(sprintf("DA %s: extra columns", b), paste(setdiff(got, want), collapse = ", "))
  } else {
    assert(sprintf("DA %s: exact column schema", b), identical(got, want),
           sprintf("missing: [%s]; unexpected: [%s]; got order: [%s]",
                   paste(setdiff(want, got), collapse = ", "),
                   paste(setdiff(got, want), collapse = ", "),
                   paste(got, collapse = ", ")))
  }
}

# ── 3. completeness ───────────────────────────────────────────────────────────
# Manifest rebuilt from OME_TISSUE_CODE exactly as qc_norm_freeze_tests.R does, plus the
# data_details string each ome's DA actually uses. Deliberate divergence from upstream, which
# validates no metab DA at all: all 48 rows are expected to produce a DA file here.
otc <- load_preflight("OME_TISSUE_CODE")

if (!is.null(otc)) {
  da_details_for <- function(ome)
    if (ome == "epigen-methylcap-seq") "malax-glmm-acute" else "dream-acute"

  man <- otc[!grepl("^lab-", otc$ome), ]
  man$token <- paste0(man$tissue_code, "_", man$ome)
  man$da_details <- vapply(man$ome, da_details_for, character(1), USE.NAMES = FALSE)
  man <- man[!duplicated(man$token), ]

  # An expected file is the token plus the right data_details: a DA written with the wrong
  # model string is a different analysis, not a rename.
  expected <- paste0(man$token, "_da_", man$da_details)
  present <- sub("_v[0-9.]+\\.txt$", "", sub("^human-precovid-sed-adu_", "", basename(da_files)))

  # No allowlist for stems that have not been rerun: this suite assumes step 09 has been run
  # to completion, so a missing DA file is a missing DA file.
  assert_completeness("every OME_TISSUE_CODE row has a DA file", present, expected)

  # file_versions.json decides what version each DA file is written as, so a token missing
  # from it gets no version and a token that is only in it is a key nobody will ever write.
  # The DA entries are the data_category == "da" leaves of its files tree — read through
  # freeze_files() rather than by grepping "_da_" out of composed filenames, so a rename in
  # any other field cannot quietly change which entries this compares.
  fv <- tryCatch({
    source(file.path(.root, "scripts", "10_build_data", "lib", "file_versions.R"))
    freeze_files()
  }, error = function(e) NULL)
  if (!is.null(fv)) {
    da <- fv[!is.na(fv$data_category) & fv$data_category == "da", ]
    fv_tokens <- paste0(da$tissue_code, "_", da$ome)
    assert_completeness("every manifest token has a file_versions.json DA key",
                        fv_tokens, man$token)
    assert_completeness("every file_versions.json DA key is a manifest token",
                        man$token, fv_tokens)
  } else skip("file_versions.json DA keys match the manifest", "config/file_versions.json unreadable")
} else skip("DA completeness", "OME_TISSUE_CODE preflight object missing")

# ── 4. per-file sanity ────────────────────────────────────────────────────────
# The grid identity nrow == n_features * n_contrasts holds by construction, not by luck:
# .convert_dream_output() loops over contrasts calling topTable(fit, coef, number = Inf,
# p.value = 1) against the SAME fit object, so every contrast returns exactly the rows the fit
# has. Measured on all 40 files in the freeze: zero ragged features.

STAT_COLS <- c("z.std", "logFC", "degrees_of_freedom", "logLik",
               "t", "AveExpr", "p_value", "adj_p_value", "methylation_diff")

for (f in da_files) {
  b <- basename(f); ome <- ome_of(f)
  if (file.size(f) / 1048576 > MAX_MB) {
    skip(sprintf("DA %s: row-level checks", b),
         sprintf("%.0f MB > DA_TESTS_MAX_MB=%s", file.size(f) / 1048576, MAX_MB))
    next
  }
  d <- fread(f, sep = "\t", header = TRUE, data.table = FALSE)
  if (!all(c("feature_id", "contrast") %in% names(d))) next   # already FAILed in section 2

  n_features <- dplyr::n_distinct(d$feature_id)
  n_contrasts <- dplyr::n_distinct(d$contrast)

  assert_unique(sprintf("DA %s: one row per feature x contrast", b),
                paste(d$feature_id, d$contrast))
  # methylcap is exempt: generate_methylcap_da.R copies Yongchao Ge's MALAX GLMM tables into
  # the freeze verbatim, so an incomplete grid there is a property of an external artifact this
  # pipeline neither fits nor can fix. Reported, not asserted — the gap stays visible without
  # failing the build for a file we only pass through.
  grid_desc <- sprintf("%d rows vs %d features x %d contrasts — %d missing",
                       nrow(d), n_features, n_contrasts, n_features * n_contrasts - nrow(d))
  if (ome == "epigen-methylcap-seq") {
    report(sprintf("DA %s: feature x contrast grid (external, not asserted)", b), grid_desc)
  } else {
    assert(sprintf("DA %s: complete feature x contrast grid", b),
           nrow(d) == n_features * n_contrasts, grid_desc)
  }

  # Contrast COUNT varies legitimately by tissue and ome (33 blood, 21 muscle/adipose, 9 for
  # adipose proteomics), so the shape is reported rather than asserted against a number.
  report(sprintf("DA %s: shape", b),
         sprintf("%d features x %d contrasts = %d rows", n_features, n_contrasts, nrow(d)))

  for (cc in intersect(STAT_COLS, names(d)))
    assert(sprintf("DA %s: %s is numeric", b, cc), is.numeric(d[[cc]]),
           sprintf("got %s", class(d[[cc]])[1]))

  in_unit_interval <- function(x) all(is.finite(x) & x >= 0 & x <= 1)
  assert(sprintf("DA %s: p_value in [0,1]", b), in_unit_interval(d$p_value),
         sprintf("%d value(s) outside [0,1] or non-finite", sum(!(is.finite(d$p_value) & d$p_value >= 0 & d$p_value <= 1))))
  assert(sprintf("DA %s: adj_p_value in [0,1]", b), in_unit_interval(d$adj_p_value),
         sprintf("%d value(s) outside [0,1] or non-finite", sum(!(is.finite(d$adj_p_value) & d$adj_p_value >= 0 & d$adj_p_value <= 1))))
  # BH can only move a p-value up
  worse <- sum(d$adj_p_value < d$p_value, na.rm = TRUE)
  assert(sprintf("DA %s: adj_p_value >= p_value", b), worse == 0,
         sprintf("%d row(s) adjusted downward", worse))

  effect_col <- if ("logFC" %in% names(d)) "logFC" else "methylation_diff"
  assert(sprintf("DA %s: %s finite", b, effect_col), all(is.finite(d[[effect_col]])),
         sprintf("%d non-finite", sum(!is.finite(d[[effect_col]]))))
  # Confidence intervals, on every dream table (methylcap is copied in verbatim and has
  # none). topTable's own CI.L/CI.R are still asserted absent, so a revert to
  # confint = TRUE fails here rather than shipping bounds built from the wrong df.
  assert(sprintf("DA %s: no uncorrected confidence-interval columns", b),
         !any(c("CI.L", "CI.R") %in% names(d)),
         paste(intersect(c("CI.L", "CI.R"), names(d)), collapse = ", "))
  ci_cols <- c("CI.L_calculated", "CI.R_calculated")
  # Gated on the columns being there, not just on this being a dream table: the checks
  # below index them directly, and a stale table from a pre-CI run would abort the whole
  # suite on an undefined-column error instead of reporting the one FAIL above.
  if ("logFC" %in% names(d)) {
    assert(sprintf("DA %s: corrected interval columns present", b),
           all(ci_cols %in% names(d)),
           paste(setdiff(ci_cols, names(d)), collapse = ", "))
  }
  if ("logFC" %in% names(d) && all(ci_cols %in% names(d))) {
    nonfinite <- sum(!is.finite(as.matrix(d[, ci_cols])))
    assert(sprintf("DA %s: interval bounds finite", b), nonfinite == 0,
           sprintf("%d non-finite bound(s)", nonfinite))
    # The bounds bracket logFC and are centred on it: both are built as logFC -/+ the same
    # non-negative margin, so a sign flip or a shifted centre means the margin came from
    # the wrong row.
    outside <- sum(!(d$CI.L_calculated <= d$logFC & d$logFC <= d$CI.R_calculated))
    assert(sprintf("DA %s: CI.L_calculated <= logFC <= CI.R_calculated", b), outside == 0,
           sprintf("%d row(s) with logFC outside its interval", outside))
    centre_gap <- max(abs((d$CI.L_calculated + d$CI.R_calculated) / 2 - d$logFC))
    assert(sprintf("DA %s: interval centred on logFC", b),
           centre_gap < 1e-6, sprintf("centre off by up to %.3g", centre_gap))
    # The margin is se * qt(0.975, df) and se is logFC/t, so dividing the half-width back
    # by se must recover a t quantile: >= qnorm(0.975) = 1.959964, which qt(0.975, df)
    # approaches from above as df grows. A margin built from a normal quantile, a
    # different alpha, or the wrong row fails this. Rows where logFC is 0 give se = 0 and
    # a 0/0 quantile, so they are excluded rather than asserted on.
    se <- d$logFC / d$t
    implied_q <- ((d$CI.R_calculated - d$CI.L_calculated) / 2) / se
    implied_q <- implied_q[is.finite(implied_q)]
    bad_q <- sum(implied_q < qnorm(0.975) - 1e-8)
    assert(sprintf("DA %s: half-width is a 95%% t quantile of the standard error", b),
           bad_q == 0,
           sprintf("%d row(s) with implied quantile below %.6f (min %.6f)",
                   bad_q, qnorm(0.975), if (length(implied_q)) min(implied_q) else NA_real_))
  }

  assert(sprintf("DA %s: single full_model", b), dplyr::n_distinct(d$full_model) == 1,
         sprintf("%d distinct formulas", dplyr::n_distinct(d$full_model)))

  # Identity is ome-dependent: metab platforms carry the family in assay and the platform in
  # its own column, every other stem puts the ome straight into assay.
  if (is_metab_platform(ome)) {
    assert(sprintf("DA %s: assay = metab, platform = %s", b, ome),
           dplyr::n_distinct(d$assay) == 1 && d$assay[1] == "metab" &&
             dplyr::n_distinct(d$platform) == 1 && d$platform[1] == ome,
           sprintf("got assay = [%s], platform = [%s]",
                   paste(unique(d$assay), collapse = ", "),
                   paste(unique(d$platform), collapse = ", ")))
  } else {
    assert(sprintf("DA %s: assay = %s", b, ome),
           dplyr::n_distinct(d$assay) == 1 && d$assay[1] == ome,
           sprintf("got [%s]", paste(unique(d$assay), collapse = ", ")))
  }
}

# ── 5. qc-norm vs DA feature accounting ───────────────────────────────────────
# The truncation guard, and the reason this file exists. Everything here WARNS or REPORTS and
# nothing fails the build.
#
# The reference is the qc-norm freeze that step 06 generated — the matrix each DA was actually
# fit over — not a pasted snapshot of released counts. A literal baseline table goes stale the
# moment either tier moves and says nothing at all about a stem nobody has rerun; the qc-norm
# freeze is regenerated alongside the DA and is always in step with it.


qc_files <- list.files(.T_FREEZE, pattern = "_qc-norm_.*\\.txt$", recursive = TRUE, full.names = TRUE)
feat_files <- list.files(.T_FREEZE, pattern = "_metadata_features_.*\\.txt$", recursive = TRUE, full.names = TRUE)
# Pair on the version-stripped token: qc-norm and DA version independently (qc-norm v1.2
# alongside DA v2.0), so the version cannot be part of the key.
qc_by_token <- stats::setNames(qc_files, sub("_qc-norm_.*$", "", sub("^human-precovid-sed-adu_", "", basename(qc_files))))
feat_by_token <- stats::setNames(feat_files, sub("_metadata_features_.*$", "", sub("^human-precovid-sed-adu_", "", basename(feat_files))))

first_col <- function(f) as.character(fread(f, select = 1, sep = "\t", header = TRUE)[[1]])
column_or_null <- function(f, nm) {
  if (!nm %in% header_of(f)) return(NULL)
  as.character(fread(f, select = nm, sep = "\t", header = TRUE)[[1]])
}
# Which filter, if any, explains a legitimate drop between the two tiers. NA means the stem
# fits every row of its qc-norm matrix, so the counts have to match exactly.
filter_explaining <- function(ome) {
  if (ome %in% c("prot-pr", "prot-ph")) return("filter_paired_n")
  if (is_metab_stem_ome(ome))           return("remove_unnamed_metab")
  return(NA_character_)
}

for (tk in names(da_by_token)) {
  f <- da_by_token[[tk]]; ome <- ome_of(f)
  desc <- sprintf("features %s: qc-norm vs DA", tk)
  if (!tk %in% names(qc_by_token)) { skip(desc, "no qc-norm freeze file for this token"); next }

  qc_ids <- first_col(qc_by_token[[tk]])
  da_ids <- unique(as.character(fread(f, select = "feature_id", sep = "\t", header = TRUE)[[1]]))

  # The id space is the qc-norm rownames plus, where the file has them, metadata_features'
  # feature_id and refmet_name. .eliminate_redundant_metab() historically renamed metab
  # rownames to refmet_name, so a DA built before this cycle's remove_redundant_metab = FALSE
  # change carries refmet names while a new one carries raw platform ids. The union is a safe
  # superset for both.
  id_space <- qc_ids
  if (tk %in% names(feat_by_token)) {
    ff <- feat_by_token[[tk]]
    id_space <- unique(c(id_space, column_or_null(ff, "feature_id"), column_or_null(ff, "refmet_name")))
  }

  unknown_ids <- setdiff(da_ids, id_space)
  if (length(unknown_ids))
    warn(sprintf("%s: DA feature_ids outside the qc-norm id space", tk),
         sprintf("%d of %d, e.g. %s", length(unknown_ids), length(da_ids),
                 paste(utils::head(unknown_ids, 3), collapse = ", ")))

  lost <- length(qc_ids) - length(da_ids)
  explained_by <- filter_explaining(ome)
  report(desc, sprintf("qc-norm %d -> DA %d (delta %+d, retained %.3f); expected filter: %s",
                       length(qc_ids), length(da_ids), -lost,
                       length(da_ids) / length(qc_ids),
                       if (is.na(explained_by)) "none — counts should match" else explained_by))

  # No filter stands between the two tiers for this ome, so every qc-norm feature should have
  # come back from the fit. A shortfall is a lost model, which is exactly how the dream/lme4
  # defect and the SOCK-worker defect both presented.
  if (is.na(explained_by) && lost > 0)
    warn(sprintf("%s: DA is short of its qc-norm matrix", tk),
         sprintf("%d of %d feature(s) have no DA row and no filter explains it — a lost fit",
                 lost, length(qc_ids)))
}

# ── 6. diff vs the reference bucket ───────────────────────────────────────────
# diff_freeze_vs_bucket() cannot be reused: it assumes a wide feature x sample matrix, and a
# DA table is long. INFO only — value differences against the release are expected and
# traced (muscle prot-ph sits near 0.997 for reasons that go back to its input files).
diff_da_vs_bucket <- function(local) {
  desc <- sprintf("diff vs bucket: %s", basename(local))
  remote_all <- .bucket_files()
  if (length(remote_all) == 0) return(skip(desc, "no bucket listing (gsutil/bucket unavailable)"))
  rel <- sub(paste0("^", .T_FREEZE, "/"), "", local)
  key <- sub("_v[0-9.]+\\.txt$", "", basename(local))
  cand <- remote_all[sub("_v[0-9.]+\\.txt$", "", basename(remote_all)) == key &
                     grepl(paste0("/", dirname(rel), "/"), remote_all)]
  if (length(cand) != 1) return(skip(desc, sprintf("%d bucket matches on %s", length(cand), .T_BUCKET)))
  dir.create(.T_RAW, showWarnings = FALSE, recursive = TRUE)
  cached <- file.path(.T_RAW, basename(cand[1]))
  if (!file.exists(cached) &&
      system(sprintf("gsutil cp '%s' '%s'", cand[1], cached), ignore.stdout = TRUE, ignore.stderr = TRUE) != 0)
    return(skip(desc, "gsutil cp failed"))

  effect <- if ("logFC" %in% header_of(local)) "logFC" else "methylation_diff"
  want <- c("feature_id", "contrast", effect)
  o <- fread(local,  select = want, sep = "\t", header = TRUE, data.table = FALSE)
  s <- fread(cached, select = want, sep = "\t", header = TRUE, data.table = FALSE)
  ko <- paste(o$feature_id, o$contrast); ks <- paste(s$feature_id, s$contrast)
  shared <- intersect(ko, ks)
  if (length(shared) < 3) return(report(desc, "fewer than 3 shared feature x contrast keys"))
  a <- o[match(shared, ko), effect]; b <- s[match(shared, ks), effect]
  ok <- is.finite(a) & is.finite(b)
  r <- if (sum(ok) >= 3) stats::cor(a[ok], b[ok]) else NA_real_
  report(sprintf("%s (n=%d shared)", desc, length(shared)),
         sprintf("%s cor = %.4f, max |diff| = %.4g; features built %d vs released %d",
                 effect, r, max(abs(a[ok] - b[ok])),
                 dplyr::n_distinct(o$feature_id), dplyr::n_distinct(s$feature_id)))
}
for (tk in names(da_by_token)) {
  f <- da_by_token[[tk]]
  if (file.size(f) / 1048576 > MAX_MB) {
    skip(sprintf("diff vs bucket: %s", basename(f)),
         sprintf("%.0f MB > DA_TESTS_MAX_MB=%s", file.size(f) / 1048576, MAX_MB))
    next
  }
  diff_da_vs_bucket(f)
}

finish()
