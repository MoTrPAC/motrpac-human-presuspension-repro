#!/usr/bin/env Rscript
# Tests for the *_SUM_STATS objects (step 11). Downstream these are the public aggregate
# layer: one row per (randomGroupCode, feature_id, Timepoint) carrying Count / Mean / SD,
# restricted to the acute bout and to features present in the step-09 DA freeze.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
suppressWarnings(suppressMessages({library(dplyr); library(data.table)}))
cat("== sum stats tests ==\n")

objs <- sub("\\.rda$", "", list.files(.T_DATA, pattern = "_SUM_STATS\\.rda$"))
assert("sum stats objects present", length(objs) > 0, "none built — run step 11")
if (length(objs) == 0) finish()

REQUIRED <- c("randomGroupCode", "feature_id", "Timepoint", "Count", "Mean", "SD",
              "tissue", "assay")

# The metabolomics stack: one {TISSUE}_METAB_SUM_STATS per tissue holding every research
# platform, which is the object key step 10 uses for the *_DA objects. metab-t-clinical is
# not part of it and is its own object, on this tier and on the DA tier alike — though it is
# labelled the same way, see below.
is_metab_stack <- function(nm) grepl("^[A-Z]+_METAB_SUM_STATS$", nm)

# Every metabolomics object, stack or clinical, carries the platform column. Which object a
# platform lands in and what its rows say they are are different questions: metab-t-clinical
# is its own object AND reads assay = "metab", which is what its *_DA object carries and what
# da_assemble_tests.R asserts for it by name.
is_metab_family <- function(nm) grepl("^[A-Z]+_METAB_", nm)

# The schema this tier is named and ordered by, which is the *_DA schema restricted to the
# columns the two share. A metabolomics object carries assay = "metab" and names the platform
# in its own column; every other object names its ome in assay and has no platform column.
# Derived from the object name, so a builder that wrote a platform into `assay` — or shipped
# a research platform as its own object again — is caught rather than described.
sum_stats_schema <- function(nm) {
  c("tissue", "assay", if (is_metab_family(nm)) "platform",
    "randomGroupCode", "Timepoint", "feature_id", "Count", "Mean", "SD")
}

# The DA is the source of truth for which features may appear, and it is the STEP-10
# assembled objects — the same source all_group_stats.R now reads, and what upstream's
# load_differential_analysis() returns. Re-derived here rather than imported from the builder
# so this test can still catch a builder mistake.
#
# It must not be re-derived from the step-09 freeze, which is what this block used to do.
# For metabolomics the freeze is in a different feature namespace: step 09 fits with
# `remove_redundant_metab = FALSE` and defers the CV-based dedup, so it carries native
# metabolite names ("Succinate") over the un-deduped feature space, while step 10 collapses
# to the RefMet/lowest-CV names ("Succinic acid") the *_QC objects use. Checking RefMet-keyed
# output against a native-keyed target failed every metab ome for the wrong reason.
#
# Epigen still comes from the freeze, because step 10 does not assemble it. The old pattern
# here (`_da_dream-`) also silently missed the methylcap tier entirely, which is fit by an
# external MALAX GLMM and lands as `_da_malax-glmm-*`; anchoring on `_da_<model>_v<x.y>.txt`
# picks up both. No adj_p_value < 0.05 filter is applied — the builder's epigen features are
# a subset of the unfiltered DA, so the subset assertion below still holds.
otc <- load_preflight("OME_TISSUE_CODE")
da_objs <- list.files(.T_DATA, pattern = "_DA\\.rda$", full.names = TRUE)
da_objs <- da_objs[!grepl("SPLICING_DA\\.rda$", da_objs)]
epigen_files <- list.files(file.path(.root, "staging", "freeze"),
                           pattern = "_da_.*_v[0-9.]+\\.txt$", recursive = TRUE, full.names = TRUE)
epigen_files <- epigen_files[grepl("_epigen-", basename(epigen_files))]

da <- if (length(da_objs) && !is.null(otc))
  dplyr::bind_rows(c(
    lapply(da_objs, function(f) {
      e <- new.env(); load(f, envir = e); o <- as.data.frame(e[[ls(e)[1]]])
      # metab objects stack platforms under assay = "metab" with the platform in its own
      # column, so the platform is the ome wherever there is one
      data.frame(feature_id = as.character(o$feature_id),
                 tissue = as.character(o$tissue),
                 ome = if ("platform" %in% colnames(o)) as.character(o$platform)
                       else as.character(o$assay),
                 stringsAsFactors = FALSE)
    }),
    lapply(epigen_files, function(f) {
      # fread, not read.csv: a COLUMN SUBSET — the epigen DA tier is ~7 GB across 5 files.
      d <- data.table::fread(f, sep = "\t", header = TRUE, select = "feature_id",
                             data.table = FALSE)
      tok <- sub("_da_.*$", "", sub("^human-precovid-sed-adu_", "", basename(f)))
      ome <- sub("^[^_]+_", "", tok)
      hit <- otc[otc$ome == ome & otc$tissue_code == sub("_.*$", "", tok), ]
      data.frame(feature_id = as.character(d$feature_id), ome = ome,
                 tissue = if (nrow(hit)) as.character(hit[["tissue"]][1]) else NA_character_,
                 stringsAsFactors = FALSE)
    }))) else NULL

for (nm in objs) {
  d <- load_built(nm)
  if (is.null(d)) { skip(sprintf("%s: loads", nm), "object missing"); next }

  # Both of these are checked on the object AS SAVED, before as.data.frame() below quietly
  # fixes them: it drops the grouping and, for a caller who reads the .rda directly, does not
  # run at all.
  #
  # A names attribute on either key column is a vial label per row — sample identifiers riding
  # inside the tier that exists because sample-level data cannot be released. It arrives from
  # subsetting a vialLabel-keyed lookup, so it is dropped in the builder at that line.
  assert(sprintf("%s: no sample identifiers in attributes", nm),
         is.null(names(d$randomGroupCode)) && is.null(names(d$Timepoint)),
         sprintf("%d name(s) on randomGroupCode, %d on Timepoint",
                 length(names(d$randomGroupCode)), length(names(d$Timepoint))))
  # A grouped object silently changes what n(), count() and summarize() return downstream.
  assert(sprintf("%s: shipped ungrouped", nm), !dplyr::is_grouped_df(d),
         sprintf("grouped by %s", paste(dplyr::group_vars(d), collapse = " x ")))

  d <- as.data.frame(d)

  assert_cols(nm, d, REQUIRED)
  if (!all(REQUIRED %in% names(d))) next

  # one row per group x feature x timepoint
  assert_unique(sprintf("%s: one row per group x feature x timepoint", nm),
                paste(d$randomGroupCode, d$feature_id, d$Timepoint))

  assert(sprintf("%s: Count >= 1", nm), all(d$Count >= 1, na.rm = TRUE),
         sprintf("min Count = %s", min(d$Count, na.rm = TRUE)))
  # Count is a number of samples, so it is whole. A fractional Count means something other
  # than n() reached the column.
  assert(sprintf("%s: Count is a whole number", nm),
         is.numeric(d$Count) && !anyNA(d$Count) && all(d$Count == round(d$Count)),
         sprintf("class %s, %d NA", class(d$Count)[1], sum(is.na(d$Count))))
  assert(sprintf("%s: Mean finite", nm), all(is.finite(d$Mean)),
         sprintf("%d non-finite", sum(!is.finite(d$Mean))))
  # SD is NA exactly when a group x timepoint has a single observation — both directions,
  # since either one alone also passes for an SD column that is entirely NA or never NA.
  bad_sd <- sum(is.na(d$SD) & d$Count > 1)
  assert(sprintf("%s: SD present whenever Count > 1", nm), bad_sd == 0,
         sprintf("%d row(s) with Count > 1 but NA SD", bad_sd))
  bad_sd1 <- sum(!is.na(d$SD) & d$Count == 1)
  assert(sprintf("%s: SD absent whenever Count == 1", nm), bad_sd1 == 0,
         sprintf("%d row(s) with Count == 1 but non-NA SD", bad_sd1))
  assert(sprintf("%s: SD non-negative", nm), all(d$SD >= 0, na.rm = TRUE))

  if (is.factor(d$Timepoint))
    assert_subset(sprintf("%s: Timepoint levels valid", nm), levels(droplevels(d$Timepoint)),
                  TIMEPOINT_LEVELS)
  assert(sprintf("%s: single tissue/assay", nm),
         dplyr::n_distinct(d$tissue) == 1 && dplyr::n_distinct(d$assay) == 1)

  # Named and ordered like the DA. Checked as a sequence, not a set: the point of the
  # column order is that the two tiers can be read side by side, which a set comparison
  # would pass while the columns were in any order at all.
  want <- sum_stats_schema(nm)
  assert(sprintf("%s: columns named and ordered like the DA", nm),
         identical(names(d), want),
         sprintf("got [%s], want [%s]", paste(names(d), collapse = ", "),
                 paste(want, collapse = ", ")))
  if ("platform" %in% names(d)) {
    assert(sprintf("%s: every row is assay metab, named by platform", nm),
           all(d$assay == "metab") && !anyNA(d$platform) &&
             all(grepl("^metab-", as.character(d$platform))),
           sprintf("assay [%s], platform [%s]",
                   paste(unique(d$assay), collapse = ", "),
                   paste(unique(d$platform), collapse = ", ")))
    # Both of these are about the STACK, not about carrying a platform column.
    # BLOOD_METAB_T_CLINICAL_SUM_STATS carries one too, and holds exactly one platform.
    if (is_metab_stack(nm)) {
      # The stack exists to hold more than one platform. An object carrying a single
      # platform is a build that fell back to one object per platform under the stack's name.
      assert(sprintf("%s: stacks >1 platform", nm), dplyr::n_distinct(d$platform) > 1,
             sprintf("only [%s]", paste(unique(d$platform), collapse = ", ")))
      # metab-t-clinical is its own object. In the stack it would duplicate the five analytes
      # it shares with the research platforms — Cortisol, Glycerol, KET, NEFA and Glucose.
      assert(sprintf("%s: no clinical chemistry in the stack", nm),
             !"metab-t-clinical" %in% as.character(d$platform))
    } else {
      assert(sprintf("%s: single platform outside the stack", nm),
             dplyr::n_distinct(d$platform) == 1,
             sprintf("[%s]", paste(unique(d$platform), collapse = ", ")))
    }
  }

  # every feature must come from the DA for this tissue x ome. The ome is the platform
  # wherever there is one, on this side exactly as on the DA side above — the two now
  # agree on that, which is the whole point of carrying the platform column here. Checked
  # per platform, not over the object: a feature measured on one platform and absent from
  # another's DA passes an object-wide subset check that it should fail.
  if (!is.null(da)) {
    omes <- if ("platform" %in% names(d)) unique(as.character(d$platform))
            else as.character(d$assay[1])
    for (ome in omes) {
      rows <- if ("platform" %in% names(d)) d[as.character(d$platform) == ome, ] else d
      keyed <- da[da$tissue == d$tissue[1] & da$ome == ome, "feature_id"]
      assert_subset(sprintf("%s: %s features come from its DA", nm, ome),
                    unique(rows$feature_id), unique(keyed))
    }
  } else skip(sprintf("%s: features come from its DA", nm), "no DA freeze found")
}

# diff each object against the package version (REPORT only — a reproduction difference is
# logged for review, never a FAIL). Three axes, because a difference in each means something
# different:
#
#   coverage   which (group x feature x timepoint) rows exist at all. Moves when the DA
#              feature set this step intersects against moves, or when an ome gains samples.
#   sample n   Count, and the group size it implies per (group x timepoint). Diffed BEFORE the
#              values: a Mean is a statistic OVER a set of samples, so one that agrees while
#              Count disagrees is a coincidence, not a reproduction, and the Mean alone cannot
#              tell those two apart.
#   values     Mean and SD over the shared rows.
max_or_na <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

# Sample number per (group x timepoint). Count varies across features within a cell only where
# a feature is NA in some of its samples, so the largest Count in a cell is the number of
# samples that cell was summarised from — the design-level sample number.
group_sizes <- function(x) {
  a <- x %>%
    dplyr::group_by(g = as.character(randomGroupCode), t = as.character(Timepoint)) %>%
    dplyr::summarize(n = max(Count), .groups = "drop")
  stats::setNames(a$n, paste(a$g, a$t))
}

# The package counterpart of a built object. For the metabolomics stack that is the object of
# the same name once the package ships one; until then it is the package's per-platform
# {TISSUE}_METAB_{PLATFORM}_SUM_STATS objects bound together. Those hold the same rows under
# the old layout, so the three diffs below keep reporting on the reproduction instead of
# reporting the reshape as a wholesale disappearance.
pkg_counterpart <- function(nm) {
  rda <- .find_rda(nm)
  if (!is.na(rda)) return(as.data.frame(.load_rda(rda, nm)))
  if (!is_metab_stack(nm)) return(NULL)
  pat <- sprintf("^%s_METAB_[A-Z_]+_SUM_STATS[.]rda$", sub("_METAB_SUM_STATS$", "", nm))
  parts <- list()
  for (p in .T_PKGS) {
    for (f in list.files(file.path(p, "data"), pattern = pat, full.names = TRUE)) {
      obj <- sub("[.]rda$", "", basename(f))
      # clinical chemistry is its own object under both layouts, never part of the stack
      if (obj == "BLOOD_METAB_T_CLINICAL_SUM_STATS") next
      if (is.null(parts[[obj]])) parts[[obj]] <- as.data.frame(.load_rda(f, obj))
    }
  }
  if (!length(parts)) return(NULL)
  dplyr::bind_rows(parts)
}

diff_one <- function(nm) {
  base <- sprintf("diff sum stats vs package: %s", nm)
  pkg <- pkg_counterpart(nm); if (is.null(pkg)) return(report(base, "no package .rda"))
  built <- as.data.frame(load_built(nm))
  if (!all(c("Count", "Mean", "SD") %in% names(pkg)))
    return(report(base, sprintf("package object carries no Count/Mean/SD — [%s]",
                                paste(names(pkg), collapse = ", "))))
  k <- function(x) paste(x$randomGroupCode, x$feature_id, x$Timepoint)
  kb <- k(built); kp <- k(pkg); sh <- intersect(kb, kp)
  gb <- group_sizes(built); gp <- group_sizes(pkg)

  report(sprintf("%s: coverage", base),
         sprintf("rows %d vs %d | features %d vs %d | cells %d vs %d | shared rows %d (built-only %d, package-only %d)",
                 nrow(built), nrow(pkg),
                 dplyr::n_distinct(built$feature_id), dplyr::n_distinct(pkg$feature_id),
                 length(gb), length(gp), length(sh),
                 length(setdiff(kb, kp)), length(setdiff(kp, kb))))

  ne <- intersect(names(gb), names(gp))
  ne <- ne[gb[ne] != gp[ne]]
  report(sprintf("%s: samples per group x timepoint", base),
         sprintf("%d shared cell(s), %d differ%s | built-only cell(s) %d, package-only %d",
                 length(intersect(names(gb), names(gp))), length(ne),
                 if (length(ne)) sprintf(" (e.g. %s)", paste(sprintf("%s %s vs %s",
                   utils::head(ne, 3), gb[utils::head(ne, 3)], gp[utils::head(ne, 3)]),
                   collapse = "; ")) else "",
                 length(setdiff(names(gb), names(gp))), length(setdiff(names(gp), names(gb)))))

  if (length(sh) < 3)
    return(report(sprintf("%s: per-row", base), "too few shared group x feature x timepoint"))
  b <- built[match(sh, kb), ]; p <- pkg[match(sh, kp), ]

  # %s, not %d, for everything derived from Count: a built object whose Count is double-typed
  # is caught by the assertion above, and this report must still print rather than abort the
  # run on an sprintf format error.
  cd <- abs(b$Count - p$Count)
  report(sprintf("%s: Count over %d shared row(s)", base, length(sh)),
         sprintf("%d row(s) differ (%.2f%%), max |diff| = %s | samples summarised %s vs %s",
                 sum(cd > 0, na.rm = TRUE), 100 * mean(cd > 0, na.rm = TRUE), max_or_na(cd),
                 sum(b$Count, na.rm = TRUE), sum(p$Count, na.rm = TRUE)))

  md <- abs(b$Mean - p$Mean)
  report(sprintf("%s: Mean", base),
         sprintf("max |diff| = %.4g, median |diff| = %.4g, %d row(s) > 1e-6",
                 max_or_na(md), stats::median(md, na.rm = TRUE),
                 sum(md > 1e-6, na.rm = TRUE)))

  sdd <- abs(b$SD - p$SD)
  report(sprintf("%s: SD", base),
         sprintf("max |diff| = %.4g, %d row(s) > 1e-6, %d NA-pattern mismatch(es)",
                 max_or_na(sdd), sum(sdd > 1e-6, na.rm = TRUE),
                 sum(is.na(b$SD) != is.na(p$SD))))
}
for (nm in objs) diff_one(nm)
finish()
