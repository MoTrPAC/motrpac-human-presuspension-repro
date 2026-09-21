#!/usr/bin/env Rscript
# Tests for CONTRAST_CONVERTER + the *_DA objects (step 10). These are the package-facing
# differential-analysis tier: step 09 fits the models and writes the freeze, step 10 only
# reshapes it. So the governing property is CONSERVATION — every freeze row this step reads
# must appear exactly once in exactly one object, with its statistics untouched. Most of what
# follows checks that, plus the contrast metadata joined on the way through.
#
# The expected object inventory is re-derived from the freeze here rather than imported from
# differential_analysis_results.R, so a mistake in the builder's grouping/naming shows up as
# a failure instead of being mirrored by the test.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
suppressWarnings(suppressMessages({library(dplyr); library(data.table)}))
cat("== da assemble tests ==\n")

# ── expected inventory, re-derived from the freeze ────────────────────────────
# Mirrors the builder's rules independently: every DA table the glob matches is assembled,
# metab platforms stack per tissue, clinical chemistry stays SPLIT into its two fitted tables
# rather than fused into one CLIN_CHEMISTRY_DA, and epigen is one object per freeze table.
otc <- load_preflight("OME_TISSUE_CODE")
da_files <- list.files(.T_FREEZE, pattern = "_da_.*_v[0-9.]+\\.txt$",
                       recursive = TRUE, full.names = TRUE)
assert("DA freeze present", length(da_files) > 0, "no *_da_*.txt under staging/freeze")
if (length(da_files) == 0 || is.null(otc)) finish()

inv <- do.call(rbind, lapply(da_files, function(f) {
  tok <- sub("_da_.*$", "", sub("^human-precovid-sed-adu_", "", basename(f)))
  code <- sub("_.*$", "", tok); ome <- sub("^[^_]+_", "", tok)
  hit <- otc[otc$ome == ome & otc$tissue_code == code, ]
  data.frame(file = f, ome = ome,
             tissue = if (nrow(hit)) as.character(hit[["tissue"]][1]) else NA_character_,
             stringsAsFactors = FALSE)
}))
assert("every DA file maps to a tissue via OME_TISSUE_CODE", !anyNA(inv$tissue),
       sprintf("unmapped: %s", paste(utils::head(basename(inv$file[is.na(inv$tissue)]), 3), collapse = ", ")))

inv <- inv[!is.na(inv$tissue), , drop = FALSE]
inv$ome_group <- ifelse(grepl("^metab-", inv$ome) & inv$ome != "metab-t-clinical", "metab", inv$ome)
inv$object <- sprintf("%s_%s_DA", toupper(inv$tissue),
                      vapply(inv$ome_group, function(o)
                        if (o == "transcript-rna-seq") "TRNSCRPT" else toupper(gsub("-", "_", o)),
                        character(1)))
expected <- sort(unique(inv$object))

built <- sort(sub("\\.rda$", "", list.files(.T_DATA, pattern = "_DA\\.rda$")))
# Other steps also write *_DA objects into the same data/ dir (SPLICING_DA is a 00_preflight
# leaf), so compare against this step's expected set rather than everything ending in _DA.
built <- intersect(built, c(expected, "CLIN_CHEMISTRY_DA"))
assert_completeness("all expected *_DA objects built", built, expected)
assert("no unexpected *_DA object from this step", length(setdiff(built, expected)) == 0,
       sprintf("extra: %s", paste(setdiff(built, expected), collapse = ", ")))

# Epigen IS assembled, one object per freeze table, and is covered by everything below.
# Asserted explicitly because it is the one ome group whose two model families carry different
# statistics, so "the glob happened to match" is not enough.
assert("epigen assembled, one object per freeze table",
       length(list.files(.T_DATA, pattern = "EPIGEN.*_DA\\.rda$")) ==
         sum(grepl("^epigen-", inv$ome)),
       sprintf("%d object(s) vs %d freeze table(s)",
               length(list.files(.T_DATA, pattern = "EPIGEN.*_DA\\.rda$")),
               sum(grepl("^epigen-", inv$ome))))

# ── CONTRAST_CONVERTER ────────────────────────────────────────────────────────
cc <- load_built("CONTRAST_CONVERTER")
if (is.null(cc)) {
  skip("CONTRAST_CONVERTER: loads", "object missing")
} else {
  cc <- as.data.frame(cc)
  assert_cols("CONTRAST_CONVERTER", cc,
              c("contrast_order", "contrast", "contrast_short", "contrast_type",
                "contrast_category", "randomGroupCode", "Timepoint"))
  assert("CONTRAST_CONVERTER: 33 contrasts", nrow(cc) == 33, sprintf("got %d", nrow(cc)))
  assert_unique("CONTRAST_CONVERTER: contrast unique", as.character(cc$contrast))
  # contrast_order is the vendored file's row order and sets the factor levels every *_DA
  # object inherits; if it is ever not 1..n the ordering downstream is silently wrong.
  assert("CONTRAST_CONVERTER: contrast_order is 1..n",
         identical(as.integer(cc$contrast_order), seq_len(nrow(cc))))
  assert_levels("CONTRAST_CONVERTER: contrast_type levels", cc$contrast_type,
                c("exercise_with_controls", "exercise_no_controls", "Endur_vs_Resist",
                  "baseline", "control_only"))
  assert_subset("CONTRAST_CONVERTER: contrast_category values",
                as.character(cc$contrast_category),
                c("EE-CON", "RE-CON", "EE-EE", "RE-RE", "EE-RE", "CON-CON"))
  assert_subset("CONTRAST_CONVERTER: Timepoint levels valid",
                levels(droplevels(cc$Timepoint)), TIMEPOINT_LEVELS)
  assert("CONTRAST_CONVERTER: no NA", !anyNA(cc))
  # every contrast in the vendored file must still be one the freeze actually fits
  assert_subset("CONTRAST_CONVERTER: randomGroupCode values",
                unique(as.character(cc$randomGroupCode)),
                c("ADUEndur", "ADUResist", "ADUControl", "ADUEndur - ADUResist"))
}

# ── per-object structure + conservation ───────────────────────────────────────
IDENTITY <- c("tissue", "assay", "full_model", "contrast", "contrast_short", "contrast_type",
              "contrast_category", "randomGroupCode", "Timepoint", "feature_id")
JOINED   <- c("contrast_short", "contrast_type", "contrast_category", "randomGroupCode", "Timepoint")

# Freeze data-row count per file, counted without parsing: these tables total ~890 MB and
# only the line count is wanted. wc -l then minus the header.
freeze_rows <- function(f) {
  n <- suppressWarnings(system2("wc", c("-l", shQuote(f)), stdout = TRUE, stderr = FALSE))
  n <- suppressWarnings(as.integer(sub("^\\s*([0-9]+).*$", "\\1", n[1])))
  if (is.na(n)) NA_integer_ else n - 1L
}

# The metab stack is NOT row-conserving, by design: step 09 fits with
# `remove_redundant_metab = FALSE` and defers the CV-based dedup downstream, so step 10
# collapses each platform onto the RefMet/lowest-CV feature space on the way in. Conservation
# is therefore asserted against the freeze rows that SURVIVE that filter, re-derived here
# straight from METABOLOMICS_CVS rather than imported from the builder.
CVS <- load_built("METABOLOMICS_CVS")
cv_keep <- function(tissue, platform) {
  if (is.null(CVS)) return(NULL)
  unique(CVS$feature_id[CVS$tissue == tissue & CVS$assay == platform & CVS$lowest_CV == "yes"])
}
# A platform every one of whose analytes is measured at a lower CV elsewhere keeps nothing and
# so contributes no rows and no platform level. Reported so an emptied platform is visible
# here as well as in the builder's log.
metab_inv <- inv[inv$ome_group == "metab", , drop = FALSE]
emptied <- if (is.null(CVS)) NULL else
  unique(metab_inv$ome[vapply(seq_len(nrow(metab_inv)), function(i)
    length(cv_keep(metab_inv$tissue[i], metab_inv$ome[i])) == 0, logical(1))])
report("metab platforms emptied by the CV dedup",
       if (is.null(CVS)) "METABOLOMICS_CVS missing"
       else if (length(emptied)) paste(emptied, collapse = ", ") else "none")
# native feature_id -> the id it should carry after the collapse (refmet_name, or itself).
# as.character before the ifelse: refmet_name is a FACTOR in METABOLOMICS_CVS, and ifelse()
# drops attributes, so a factor arm returns integer codes rather than names.
cv_map <- function(tissue, platform) {
  if (is.null(CVS)) return(NULL)
  lk <- unique(CVS[CVS$tissue == tissue & CVS$assay == platform & CVS$lowest_CV == "yes",
                   c("feature_id", "refmet_name")])
  rn <- as.character(lk$refmet_name); fid <- as.character(lk$feature_id)
  stats::setNames(ifelse(!is.na(rn) & nzchar(rn), rn, fid), fid)
}
# freeze rows for one metab platform file, after the lowest_CV filter
freeze_rows_metab <- function(f, tissue, platform) {
  keep <- cv_keep(tissue, platform); if (is.null(keep)) return(NA_integer_)
  d <- data.table::fread(f, sep = "\t", header = TRUE, select = "feature_id", data.table = FALSE)
  sum(d$feature_id %in% keep)
}

for (nm in expected) {
  d <- load_built(nm)
  if (is.null(d)) { skip(sprintf("%s: loads", nm), "object missing"); next }
  assert(sprintf("%s: is data.table", nm), data.table::is.data.table(d))
  dd <- as.data.frame(d)
  rows <- inv[inv$object == nm, , drop = FALSE]
  is_metab <- unique(rows$ome_group) == "metab"
  # metab-t-clinical stands outside the metab stack but its freeze file carries the same
  # assay = "metab" / platform = ome pair, so the column reaches its object too.
  carries_platform <- is_metab || "metab-t-clinical" %in% rows$ome

  assert_cols(nm, dd, c(IDENTITY, "p_value", "adj_p_value"))
  if (!all(IDENTITY %in% names(dd))) next
  if (carries_platform) assert(sprintf("%s: has platform column", nm), "platform" %in% names(dd))

  # column order: the identity block leads, in this exact order, with platform after assay
  lead <- c("tissue", "assay", if (carries_platform) "platform", "full_model", "contrast",
            "contrast_short", "contrast_type", "contrast_category", "randomGroupCode",
            "Timepoint", "feature_id")
  assert(sprintf("%s: leading column order", nm), identical(names(dd)[seq_along(lead)], lead),
         sprintf("got %s", paste(utils::head(names(dd), length(lead)), collapse = " ")))
  if ("z.std" %in% names(dd))
    assert(sprintf("%s: z.std sits before p_value", nm),
           which(names(dd) == "z.std") == which(names(dd) == "p_value") - 1L)

  # The two epigen model families carry different statistics and are the only place in this
  # step where that is true, so each is asserted rather than left to intersect(). methylcap is
  # the MALAX GLMM table step 09 passes through verbatim; ATAC is dream-fit like every other ome.
  ome_i <- unique(rows$ome)
  if (identical(ome_i, "epigen-methylcap-seq")) {
    assert_cols(sprintf("%s: malax statistics", nm), dd,
                c("methylation_diff", "t", "AveExpr", "p_value", "adj_p_value"))
    assert(sprintf("%s: methylation_diff, not logFC", nm),
           "methylation_diff" %in% names(dd) && !"logFC" %in% names(dd))
    # CI.L/CI.R are no longer discriminating — no table has them since run_dream() moved to
    # confint = FALSE. The corrected CI.L_calculated/CI.R_calculated pair IS dream-only,
    # but the set below is left as z.std / degrees_of_freedom / logLik: those three come
    # straight off the fit, while the interval is derived, so they fail sooner and clearer.
    assert(sprintf("%s: dream-only columns absent", nm),
           !any(c("z.std", "degrees_of_freedom", "logLik") %in% names(dd)),
           paste(intersect(c("z.std", "degrees_of_freedom", "logLik"),
                           names(dd)), collapse = ", "))
  } else if (identical(ome_i, "epigen-atac-seq")) {
    assert_cols(sprintf("%s: dream statistics", nm), dd,
                c("logFC", "z.std", "degrees_of_freedom", "logLik",
                  "t", "AveExpr", "p_value", "adj_p_value"))
    assert(sprintf("%s: logFC, not methylation_diff", nm),
           "logFC" %in% names(dd) && !"methylation_diff" %in% names(dd))
  }

  assert(sprintf("%s: single tissue", nm), dplyr::n_distinct(dd$tissue) == 1,
         paste(unique(dd$tissue), collapse = ", "))
  assert(sprintf("%s: single assay", nm), dplyr::n_distinct(dd$assay) == 1,
         paste(unique(dd$assay), collapse = ", "))
  assert(sprintf("%s: tissue matches its freeze files", nm),
         identical(unique(dd$tissue), unique(rows$tissue)))

  # CONSERVATION: one row per feature x contrast (x platform), and the row count equals the
  # sum of the freeze tables that fed it. Together these say assembly neither dropped,
  # duplicated, nor invented a row.
  key <- if (is_metab) paste(dd$platform, dd$feature_id, dd$contrast) else paste(dd$feature_id, dd$contrast)
  assert_unique(sprintf("%s: one row per %sfeature x contrast", nm, if (is_metab) "platform x " else ""), key)
  want <- if (is_metab)
    sum(vapply(seq_len(nrow(rows)),
               function(i) freeze_rows_metab(rows$file[i], unique(rows$tissue), rows$ome[i]),
               integer(1)))
  else sum(vapply(rows$file, freeze_rows, integer(1)))
  if (is.na(want)) skip(sprintf("%s: row count matches freeze", nm), "could not count freeze rows")
  else assert(sprintf("%s: row count matches freeze%s", nm, if (is_metab) " (after CV collapse)" else ""),
              nrow(dd) == want,
              sprintf("%d built vs %d in %d freeze file(s)", nrow(dd), want, nrow(rows)))

  # The metab namespace must be RefMet, not native. Catches a missed or partial collapse:
  # any feature carrying a native id that METABOLOMICS_CVS maps to a refmet_name is one the
  # rename should have replaced, and would sit in a different namespace from the *_QC objects
  # that step 11 intersects against.
  if (is_metab && !is.null(CVS)) {
    # Exact set equality per platform, against the mapping re-derived from METABOLOMICS_CVS.
    # Deliberately not a weaker "no native ids survived" check: the first version of this
    # collapse renamed every metabolite to a factor's integer codes ("79", "182"), which
    # contains no native ids and so passed that check while being entirely wrong.
    wrong <- character(0)
    for (i in seq_len(nrow(rows))) {
      mp <- cv_map(unique(rows$tissue), rows$ome[i]); if (is.null(mp)) next
      raw <- data.table::fread(rows$file[i], sep = "\t", header = TRUE,
                               select = "feature_id", data.table = FALSE)
      want <- unique(unname(mp[intersect(unique(raw$feature_id), names(mp))]))
      got  <- unique(dd$feature_id[as.character(dd$platform) == rows$ome[i]])
      if (!setequal(got, want))
        wrong <- c(wrong, sprintf("%s (built %d vs expected %d; e.g. built '%s' expected '%s')",
                                  rows$ome[i], length(got), length(want),
                                  paste(utils::head(setdiff(got, want), 2), collapse = ","),
                                  paste(utils::head(setdiff(want, got), 2), collapse = ",")))
    }
    assert(sprintf("%s: features match the RefMet/lowest-CV mapping exactly", nm),
           length(wrong) == 0,
           paste(utils::head(wrong, 3), collapse = " | "))
  }

  # the contrast join must be total — an NA here means a contrast reached the object that
  # contrast_converter.txt does not describe
  for (cl in JOINED)
    assert(sprintf("%s: %s fully joined", nm, cl), !anyNA(dd[[cl]]),
           sprintf("%d NA", sum(is.na(dd[[cl]]))))
  if (!is.null(cc))
    assert_subset(sprintf("%s: contrast levels come from CONTRAST_CONVERTER", nm),
                  levels(droplevels(dd$contrast)), as.character(cc$contrast))

  assert(sprintf("%s: feature_id is character", nm), is.character(dd$feature_id))
  assert(sprintf("%s: feature_id has no NA", nm), !anyNA(dd$feature_id))
  assert(sprintf("%s: p_value in [0,1]", nm),
         all(dd$p_value >= 0 & dd$p_value <= 1, na.rm = TRUE))
  assert(sprintf("%s: adj_p_value in [0,1]", nm),
         all(dd$adj_p_value >= 0 & dd$adj_p_value <= 1, na.rm = TRUE))
  # BH adjustment is monotone and never shrinks a p-value; a violation means the freeze was
  # adjusted across the wrong grouping or the columns were transposed on the way in.
  bad <- sum(dd$adj_p_value < dd$p_value - 1e-9, na.rm = TRUE)
  assert(sprintf("%s: adj_p_value >= p_value", nm), bad == 0, sprintf("%d row(s) below", bad))

  if (is_metab) {
    # The stack must carry exactly the platforms whose files fed it, named by their ome —
    # less any platform the CV dedup emptied, which contributes no row and so no level. That
    # set is re-derived from METABOLOMICS_CVS rather than taken from the builder.
    want_pf <- rows$ome
    if (!is.null(CVS))
      want_pf <- rows$ome[vapply(seq_len(nrow(rows)),
                                 function(i) length(cv_keep(unique(rows$tissue), rows$ome[i])) > 0,
                                 logical(1))]
    assert(sprintf("%s: platforms match its freeze files", nm),
           setequal(unique(as.character(dd$platform)), want_pf),
           sprintf("built [%s] vs expected [%s]",
                   paste(sort(unique(as.character(dd$platform))), collapse = ", "),
                   paste(sort(want_pf), collapse = ", ")))
    assert(sprintf("%s: assay is metab", nm), all(dd$assay == "metab"))
  }
}

# ── statistics survive assembly unchanged ─────────────────────────────────────
# The point of the step: it reshapes, it does not refit. So every statistic in an assembled
# object must be bit-identical to the freeze cell it came from. This is what separates a real
# assembly bug from the expected "DIFFER (content)" against the released package below —
# that difference is step 09 refitting, and it must not be reachable from here.
#
# Checked by re-reading each freeze table and joining on the object's own key, one file at a
# time so peak memory stays at one table.
STAT_COLS <- c("logFC", "t", "AveExpr", "z.std", "p_value", "adj_p_value",
               "degrees_of_freedom", "logLik", "methylation_diff")
for (nm in expected) {
  d <- load_built(nm); if (is.null(d)) next
  rows <- inv[inv$object == nm, , drop = FALSE]
  is_metab <- unique(rows$ome_group) == "metab"
  bad <- character(0); checked <- 0L
  for (i in seq_len(nrow(rows))) {
    raw <- data.table::fread(rows$file[i], sep = "\t", header = TRUE, data.table = FALSE)
    # For metab, put the freeze rows through the same CV filter + RefMet rename the builder
    # applied, so the surviving rows can be matched in the object's namespace. Rows the
    # filter drops are not expected to be there and are not checked.
    if (is_metab) {
      mp <- cv_map(unique(rows$tissue), rows$ome[i])
      if (is.null(mp)) { bad <- c(bad, sprintf("%s: no CV map", rows$ome[i])); next }
      raw <- raw[raw$feature_id %in% names(mp), , drop = FALSE]
      raw$feature_id <- unname(mp[raw$feature_id])
    }
    # A platform whose every feature loses the lowest-CV contest contributes no rows to the
    # object, and none are expected. Skip before keying: paste() recycles its scalar against
    # a zero-row frame and would yield one key matching nothing.
    if (nrow(raw) == 0L) { checked <- checked + 1L; next }
    stats <- intersect(STAT_COLS, intersect(colnames(raw), colnames(d)))
    kb <- if (is_metab) paste(as.character(d$platform), d$feature_id, as.character(d$contrast))
          else paste(d$feature_id, as.character(d$contrast))
    kr <- if (is_metab) paste(rows$ome[i], raw$feature_id, raw$contrast)
          else paste(raw$feature_id, raw$contrast)
    m <- match(kr, kb)
    if (anyNA(m)) { bad <- c(bad, sprintf("%s: %d row(s) absent from object", rows$ome[i], sum(is.na(m)))); next }
    for (s in stats) if (!identical(as.numeric(raw[[s]]), as.numeric(d[[s]][m])))
      bad <- c(bad, sprintf("%s/%s", rows$ome[i], s))
    checked <- checked + 1L
  }
  assert(sprintf("%s: statistics identical to the freeze", nm), length(bad) == 0,
         sprintf("%d mismatch(es): %s", length(bad), paste(utils::head(bad, 4), collapse = ", ")))
  rm(d); invisible(gc(verbose = FALSE))
}

# ── clinical chemistry stays split ────────────────────────────────────────────
# Deliberate divergence from the released package, which fuses these two fitted tables into
# one CLIN_CHEMISTRY_DA under the synthetic assay "clinical-chemistry". docs/data_objects.tsv
# lists them separately, so the split is asserted, and the union is reported against the
# package object as the check that nothing went missing in splitting it.
pc <- load_built("BLOOD_PROT_CLINICAL_DA"); mc <- load_built("BLOOD_METAB_T_CLINICAL_DA")
if (is.null(pc) || is.null(mc)) {
  skip("clinical chemistry split into two objects", "one or both objects missing")
} else {
  pc <- as.data.frame(pc); mc <- as.data.frame(mc)
  assert("clinical chemistry: prot-clinical and metab-t-clinical are separate objects",
         nrow(pc) > 0 && nrow(mc) > 0)
  assert("clinical chemistry: the two objects share no feature",
         length(intersect(pc$feature_id, mc$feature_id)) == 0,
         sprintf("shared: %s", paste(intersect(pc$feature_id, mc$feature_id), collapse = ", ")))
  # prot-clinical names itself in assay; metab-t-clinical follows the metab convention and
  # names itself in platform, under assay = "metab".
  assert("clinical chemistry: assays stay distinct",
         identical(unique(as.character(pc$assay)), "prot-clinical") &&
         identical(unique(as.character(mc$assay)), "metab") &&
         identical(unique(as.character(mc$platform)), "metab-t-clinical"))
  ref <- .find_rda("CLIN_CHEMISTRY_DA")
  if (is.na(ref)) report("clinical chemistry: union vs package CLIN_CHEMISTRY_DA", "no package .rda")
  else {
    cch <- as.data.frame(.load_rda(ref, "CLIN_CHEMISTRY_DA"))
    report("clinical chemistry: union vs package CLIN_CHEMISTRY_DA",
           sprintf("%d + %d = %d rows / %d features vs %d rows / %d features",
                   nrow(pc), nrow(mc), nrow(pc) + nrow(mc),
                   length(unique(c(pc$feature_id, mc$feature_id))),
                   nrow(cch), length(unique(as.character(cch$feature_id)))))
  }
}

# ── diff vs the shipped package objects (INFO only) ───────────────────────────
# Never a FAIL: step 09 refits every dream ome locally, so row counts legitimately differ
# from the released package where the freeze content has moved on.
for (nm in c("CONTRAST_CONVERTER", expected)) {
  b <- load_built(nm); rda <- .find_rda(nm)
  if (is.null(b)) next
  if (is.na(rda)) { report(sprintf("diff vs package: %s", nm), "no package .rda"); next }
  p <- .load_rda(rda, nm)
  report(sprintf("diff vs package: %s", nm),
         sprintf("%s | built %d x %d vs pkg %d x %d | cols %s",
                 .cmp(as.data.frame(b), as.data.frame(p)),
                 nrow(b), ncol(b), nrow(p), ncol(p),
                 if (identical(colnames(b), colnames(p))) "identical" else "DIFFER"))
}

finish()
