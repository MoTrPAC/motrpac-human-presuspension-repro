#!/usr/bin/env Rscript
# Stage 1 data objects: the *_SUM_STATS objects  (deps: QC_NORM, DA)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/all_group_stats.R.
#
# Sample-level expression cannot be released publicly, so this is the aggregate layer:
# per (randomGroupCode, feature_id, Timepoint) Count / Mean / SD, restricted to the acute
# bout (visitcode == "ADU_BAS") and to the features that actually appear in the differential
# analysis for that tissue x ome.
#
# Output: scripts/10_build_data/data/{TISSUE}_{OME}_SUM_STATS.rda, with the columns
#   tissue, assay, [platform], randomGroupCode, Timepoint, feature_id, Count, Mean, SD
# named, ordered AND KEYED to agree with the *_DA objects: the whole metabolomics family
# reads assay = "metab" with the platform in its own column, and every other ome names
# itself in assay and has no platform column. The research platforms stack into one
# {TISSUE}_METAB_SUM_STATS per tissue; metab-t-clinical is labelled the same way but stays
# its own object, exactly as it does on the DA tier.
#
# The stack is what step 10 does for the DA (see its ome_group), and this tier used to ship
# one object per metabolomics platform instead: fourteen {TISSUE}_METAB_{PLATFORM}_SUM_STATS
# objects against a single {TISSUE}_METAB_DA. A caller loading both then had a list nested
# one way on one side and another way on the other, and had to know which vocabulary each
# tier used before it could join them.
#
# ---- Adaptations from upstream ----------------------------------------------------------
#
# 1. DA comes from the step-10 assembled *_DA objects, which is what upstream reads too —
#    load_differential_analysis() returns exactly those. Only epigen still comes from the
#    step-09 freeze, because step 10 does not assemble it (upstream ships no epigen .rda).
#
#    This step used to read the freeze for everything, as a stand-in while step 10 was a
#    stub. That was not equivalent: for metabolomics the freeze is in a DIFFERENT FEATURE
#    NAMESPACE from the *_QC objects summarised here, so the intersection below came up empty
#    for 9 of 47 tissue x ome combos. See the DA block for the full account.
#
#    The metab objects stack every platform under assay = "metab" with the platform in its
#    own column, so the assay column does not name the ome; the platform column is used as
#    the ome wherever there is one. This step summarises one platform at a time, so it needs
#    the ome per row on the DA side however the stack is labelled.
#
# 2. Only feature_id and adj_p_value are consumed, which is also what lets the epigen tier be
#    read alongside the rest without reconciling schemas: methylcap arrives from an external
#    MALAX GLMM pipeline with 9 columns (methylation_diff as the effect size and no z.std /
#    degrees_of_freedom / logLik) against 12 for a dream table. Both families
#    carry the two columns needed here. See read_epigen_da() below.
#


suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(data.table)}))

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

out_dir <- file.path(ROOT, "scripts", "10_build_data", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- differential analysis ----------------------------------------------------------------
# Sourced from the STEP-10 assembled *_DA objects, which is what upstream does
# (load_differential_analysis() reads exactly these). This file previously read the step-09
# freeze directly, because step 10 was still a stub; that workaround is what the note here
# used to describe, and it is now wrong to keep.
#
# EPIGEN arrives through the same glob as everything else: step 10 assembles it into
# {TISSUE}_{OME}_DA.rda alongside the other objects. It must NOT also be read from the freeze
# here — that would put every epigen feature in da_data twice. It is still summarised
# restricted to adj_p_value < 0.05, below.
da_dir <- file.path(.STAGE1_DATA)
da_objs <- list.files(da_dir, pattern = "_DA\\.rda$", full.names = TRUE)
da_objs <- da_objs[!grepl("SPLICING_DA\\.rda$", da_objs)]   # a 00_preflight leaf, not this tier
if (length(da_objs) == 0)
  stop("no assembled *_DA objects in ", da_dir, " — run step 10 first")

read_da_obj <- function(f) {
  e <- new.env(); load(f, envir = e); o <- as.data.frame(e[[ls(e)[1]]])
  # The metab objects stack every platform under assay = "metab" and carry the platform in
  # its own column. This step matches per platform (curr_assay is e.g. "metab-u-rppos"), so
  # the platform is the ome wherever there is one.
  ome <- if ("platform" %in% colnames(o)) as.character(o$platform) else as.character(o$assay)
  data.frame(feature_id = as.character(o$feature_id), adj_p_value = o$adj_p_value,
             tissue = as.character(o$tissue), ome = ome, stringsAsFactors = FALSE)
}

da_data <- dplyr::bind_rows(lapply(da_objs, read_da_obj))
message(sprintf("DA: %d assembled object(s), %d rows, %d tissue x ome combo(s)",
                length(da_objs), nrow(da_data),
                nrow(unique(da_data[, c("tissue", "ome")]))))

# --- qc-norm: every ome, including clinical and epigen -------------------------------------
# Epigen is included (upstream passes epigen = FALSE). Step 08 now builds the ATAC and
# methylcap *_QC objects, so they arrive through the ordinary load. The adj_p_value < 0.05
# restriction applied to epigen below is upstream's and is kept: those matrices carry
# hundreds of thousands of features, and only the significant ones are summarised.
qc_data <- load_qc_local(selected_tissues = "all",
                         selected_omes = "all",
                         load_acute_only = FALSE,
                         load_clinical = TRUE)

# Two different questions, and they have different answers for metab-t-clinical.
#
# WHICH OBJECT a summary statistic goes in — step 10's ome_group. Every research
# metabolomics platform belongs to one "metab" object per tissue; every other ome is its own
# object. metab-t-clinical is excluded: it is clinical chemistry, published on its own, and
# putting it in the stack would put the five analytes it shares with the research platforms
# (Cortisol, Glycerol, KET, NEFA, Glucose) into one object twice.
in_metab_stack <- function(ome) grepl("^metab-", ome) & ome != "metab-t-clinical"

# WHAT THE ROWS SAY THEY ARE — the whole metabolomics family, clinical included, reads
# assay = "metab" and names its platform in the platform column. That is the freeze's own
# labelling, it is what every *_DA object carries, and da_assemble_tests.R asserts it for
# BLOOD_METAB_T_CLINICAL_DA by name. Being its own object and being labelled "metab" are
# not in tension: `platform` is what tells the two apart, on both tiers, and any key that
# has to separate clinical chemistry from a research platform includes it.
in_metab_family <- function(ome) grepl("^metab-", ome)

write_sum_stats <- function(nm, x) {
  assign(nm, x)
  save(list = nm, file = file.path(out_dir, paste0(nm, ".rda")), compress = TRUE, version = 3)
  nm
}

written <- character(0)
for (curr_tissue in names(qc_data)) {
  # The stack for this tissue, one element per platform, bound and written once the ome loop
  # below has finished with the tissue.
  metab_parts <- list()

  for (curr_assay in names(qc_data[[curr_tissue]])) {
    curr_data <- qc_data[[curr_tissue]][[curr_assay]][["qc_norm"]]
    if (is.null(curr_data) || nrow(curr_data) == 0) next

    # keep only features that appear in this tissue x assay's differential analysis
    matching_da <- da_data %>%
      dplyr::filter(.data$ome == curr_assay, .data$tissue == curr_tissue)
    if (grepl("^epigen-", curr_assay))
      matching_da <- matching_da %>% dplyr::filter(.data$adj_p_value < 0.05)
    matching_da <- unique(matching_da$feature_id)
    if (length(matching_da) == 0) {
      message(sprintf("  %-8s %-22s no DA rows — skipped", curr_tissue, curr_assay))
      next
    }
    curr_data <- curr_data[rownames(curr_data) %in% matching_da, , drop = FALSE]

    curr_meta <- qc_data[[curr_tissue]][[curr_assay]][["sample_metadata"]] %>%
      dplyr::filter(visitcode == "ADU_BAS")
    curr_data <- curr_data[, colnames(curr_data) %in% curr_meta$vialLabel, drop = FALSE]
    curr_meta <- curr_meta[match(colnames(curr_data), curr_meta$vialLabel), ]

    # Nothing survived one of the two intersections above. Reported and skipped rather than
    # left to fail: t(curr_data) on an empty frame yields a single Sample column, and the
    # pivot below then dies with "`cols` must select at least one column" — an error that
    # names neither the tissue, the ome, nor which intersection came up empty. This is how
    # the metab namespace mismatch first surfaced.
    if (nrow(curr_data) == 0 || ncol(curr_data) == 0) {
      message(sprintf("  %-8s %-22s SKIPPED — %d feature(s) x %d acute sample(s) after intersecting with DA",
                      curr_tissue, curr_assay, nrow(curr_data), ncol(curr_data)))
      next
    }

    vial_to_group     <- setNames(curr_meta$randomGroupCode, curr_meta$vialLabel)
    vial_to_timepoint <- setNames(curr_meta$Timepoint, curr_meta$vialLabel)

    curr_data_long <- as.data.frame(t(curr_data))
    curr_data_long$Sample <- rownames(curr_data_long)
    curr_data_long <- tidyr::pivot_longer(curr_data_long, -Sample,
                                          names_to = "feature_id", values_to = "Value") %>%
      dplyr::filter(!is.na(Value))
    # unname(): both lookups are keyed by vialLabel, so subsetting them carries the vial label
    # of every row along as a names attribute on the column. These objects are the public
    # aggregate layer and must not carry sample identifiers — not in a column, and not in an
    # attribute nobody prints. Whether the attribute would survive the summarize() below is a
    # dplyr implementation detail; it is dropped here, where it is created.
    curr_data_long$randomGroupCode <- unname(vial_to_group[curr_data_long$Sample])
    curr_data_long$Timepoint <- unname(vial_to_timepoint[curr_data_long$Sample])

    # .groups = "drop": summarize() defaults to "drop_last", which would ship an object still
    # grouped by (randomGroupCode, feature_id). Downstream that changes what n(), count() and
    # a further summarize() return — per group rather than over the table.
    group_stats <- curr_data_long %>%
      dplyr::group_by(randomGroupCode, feature_id, Timepoint) %>%
      dplyr::summarize(Count = dplyr::n(), Mean = mean(Value), SD = sd(Value), .groups = "drop")
    # Named the way the *_DA objects name themselves. DA stacks the research metabolomics
    # platforms under assay = "metab" and carries the platform in its own column, while
    # every other ome names itself in assay and has no platform column at all. This tier
    # used to put the platform in `assay` instead, so a summary statistic and the
    # differential analysis row it belongs to disagreed about what the ome was called,
    # and every consumer joining the two had to translate between the two vocabularies.
    # read_da_obj() above is the translation this removes the need for.
    is_metab <- in_metab_family(curr_assay)
    is_stack <- in_metab_stack(curr_assay)
    group_stats$tissue <- curr_tissue
    group_stats$assay <- if (is_metab) "metab" else curr_assay
    if (is_metab) group_stats$platform <- curr_assay

    # and ordered the way DA orders the columns it shares: what the rows are
    # (tissue, assay, platform), then what identifies a row (group, timepoint, feature),
    # then the statistics. DA puts its contrast columns between the two, and this tier
    # has none, so the shared columns line up.
    group_stats <- group_stats[, c("tissue", "assay",
                                   if (is_metab) "platform",
                                   "randomGroupCode", "Timepoint", "feature_id",
                                   "Count", "Mean", "SD"), drop = FALSE]

    # A research metabolomics platform is one part of its tissue's stack, not an object.
    if (is_stack) {
      metab_parts[[curr_assay]] <- group_stats
      message(sprintf("  %-8s %-22s %5d feature(s) -> stacked (%d rows)",
                      curr_tissue, curr_assay, length(unique(group_stats$feature_id)),
                      nrow(group_stats)))
      next
    }

    nm <- gsub("-", "_", paste(toupper(curr_tissue), toupper(curr_assay), "SUM_STATS", sep = "_"))
    written <- c(written, write_sum_stats(nm, group_stats))
    message(sprintf("  %-8s %-22s %5d feature(s) -> %s (%d rows)",
                    curr_tissue, curr_assay, length(unique(group_stats$feature_id)),
                    nm, nrow(group_stats)))
  }

  # --- the metabolomics stack for this tissue ---------------------------------------------
  if (length(metab_parts)) {
    metab_parts <- metab_parts[sort(names(metab_parts))]
    stacked <- as.data.frame(dplyr::bind_rows(metab_parts))
    nm <- paste0(toupper(curr_tissue), "_METAB_SUM_STATS")
    written <- c(written, write_sum_stats(nm, stacked))
    message(sprintf("  %-8s %-22s %5d feature(s) -> %s (%d rows, %d platform(s))",
                    curr_tissue, "metab", length(unique(stacked$feature_id)),
                    nm, nrow(stacked), length(metab_parts)))
  }
}

message(sprintf("SUM_STATS: wrote %d object(s) -> %s", length(written), out_dir))
