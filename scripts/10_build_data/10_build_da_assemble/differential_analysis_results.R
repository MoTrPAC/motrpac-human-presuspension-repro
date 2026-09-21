#!/usr/bin/env Rscript
# Stage 1 data objects: CONTRAST_CONVERTER + the *_DA objects  (deps: DA)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/differential_analysis_results.R.
#
# This is the assembly tier. Step 09 FITS the models and writes one BIC-named table per
# tissue x ome into staging/freeze/<ome-group>/da/; this step reads those tables back,
# attaches the contrast metadata, and reshapes them into the package-facing `*_DA` objects.
# Nothing is refit here — no row's statistics change, only which object it lands in, what
# the columns are called, and what order they sit in.
#
# Output: scripts/10_build_data/data/{CONTRAST_CONVERTER,{TISSUE}_{OME}_DA}.rda
#
# ---- Adaptations from upstream ----------------------------------------------------------
#
# 1. DA comes from this stage's step-09 freeze rather than
#    MotrpacHumanPreSuspensionAnalysis:::.load_differential_analysis(), which gsutil-copies
#    the released DA repo into a tempdir. The freeze IS that repo's content, already local,
#    so the download stem is dropped entirely. Everything downstream of the load — the
#    metabolomics stacking, .process_raw_DA(), the naming — is ported as-is.
#
# 2. The freeze carries `assay` but NOT `tissue` (.convert_dream_output() sets a tissue
#    column and drops it again in its closing select()). Upstream recovers both from the
#    "<tissue>.<ome>" list names its loader builds; here they are recovered from the BIC
#    tissue code in each filename via OME_TISSUE_CODE, the same way steps 08 and 11 read the
#    freeze. Note two codes map to muscle (t06-muscle, t10-muscle) and four to blood
#    (t02-plasma, t03-edta, t04-blood-rna, t05-pbmc), so the code is not the tissue.
#
# 3. OME_TISSUE_CODE is read straight from the preflight .rds instead of via lib/qc_helpers.R.
#    qc_helpers loads pheno.rda and the QC objects on source, none of which this step needs —
#    it never touches sample-level data. Reading the one lookup it does use keeps the
#    dependency honest: this step needs step 09's freeze and preflight's code table, nothing
#    else.
#
# 4. `usethis::use_data()` is replaced by a plain save() into scripts/10_build_data/data/,
#    matching every other adapted step. use_data() writes into an installed package source
#    tree; this pipeline stages objects and promotes them separately.
#
# ---- What is and is not assembled -------------------------------------------------------
#
# EPIGEN IS ASSEMBLED, unlike upstream's `epigen = FALSE`. The 5 ATAC and methylcap tables
# reshape through the same path as everything else and land beside the other objects as
# {TISSUE}_{OME}_DA.rda. Moving them to AWS, where load_differential_analysis(epigen = TRUE)
# fetches from, is handled separately.
#
# This is where the methylcap schema caveat in ../09_build_da/sources/README.md ("whoever
# ports it must reconcile the two schemas rather than assume every ome has a logFC") gets
# exercised. process_da() orders and keys columns through intersect() and never names a
# statistic, so the MALAX tables — methylation_diff instead of logFC, and none of dream's
# z.std / degrees_of_freedom / logLik — reshape rather than error.
#
# Two consumers glob `_DA\.rda$` over the output dir and both had to be told about epigen:
# step 11 now reads it from these objects instead of re-reading the freeze (it would
# otherwise count every epigen feature twice), and step 12 skips it at the glob rather than
# loading ~7 GB it drops when it filters to its five enrichment omes.
#
# CLINICAL CHEMISTRY IS SPLIT, and this is the one place the local objects deliberately
# differ from the released package. Upstream ships a single CLIN_CHEMISTRY_DA holding
# 9 features x 33 contrasts; those 9 features arrive in the freeze as TWO separately-fit
# tables — prot-clinical (CK, Glucagon, Insulin) and metab-t-clinical (Cortisol, Glucose,
# Glycerol, KET, Lactate, NEFA) — which upstream concatenates under the synthetic assay name
# "clinical-chemistry". They are kept apart here as BLOOD_PROT_CLINICAL_DA and
# BLOOD_METAB_T_CLINICAL_DA, one object per fitted table, which is what
# docs/data_objects.tsv lists. Merging them would fuse two different assays under a name
# neither one uses; keeping them split means every object in this step maps back to exactly
# one freeze file (or, for metab, one stack of platform files).
#
# metab-t-imm-crt reaches the metab stack and comes out of it empty. Its single feature is
# Cortisol — the same analyte metab-t-clinical carries, measured by immunoassay — and blood
# metab-u-hilicpos measures it at a lower CV, so METABOLOMICS_CVS flags no imm-crt row
# lowest_CV == "yes" and none survives collapse_metab_cv(). The dedup that keeps every other
# duplicated analyte to one copy is what keeps a second Cortisol row out of BLOOD_METAB_DA,
# which is the outcome upstream reaches by excluding the platform outright. An emptied
# platform is reported in the log rather than treated as an error.

suppressWarnings(suppressMessages({library(dplyr); library(data.table)}))

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
HERE    <- file.path(ROOT, "scripts", "10_build_data", "10_build_da_assemble")
out_dir <- file.path(ROOT, "scripts", "10_build_data", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- CONTRAST_CONVERTER -----------------------------------------------------------------
# Ported verbatim from upstream. The vendored file supplies only `contrast` (the model-matrix
# expression the DA tables key on) and `contrast_short`; every other column is derived here,
# so this block is the definition of contrast_type / contrast_category / randomGroupCode /
# Timepoint for the whole package. ROW ORDER IS SIGNIFICANT: it sets contrast_order and the
# factor level order that the *_DA objects inherit, which is why the file is vendored rather
# than rebuilt from the contrasts present in the freeze.
cc_file <- file.path(HERE, "sources", "contrast_converter.txt")
if (!file.exists(cc_file)) stop("missing vendored contrast converter: ", cc_file)

CONTRAST_CONVERTER <- read.csv(cc_file, sep = "\t", check.names = FALSE,
                               stringsAsFactors = FALSE) %>%
  mutate(
    contrast_type = case_when(
      grepl("^[^.]+\\.pre_exercise - [^.]+\\.pre_exercise$",
            contrast_short) ~ "baseline",
      grepl("^Control", contrast_short) ~ "control_only",
      grepl("Endur", contrast_short) &
        grepl("Resist", contrast_short) ~ "Endur_vs_Resist",
      grepl("^[^. ]+.[^ ]+ - [^.]+\\.[^ ]+$",
            contrast_short) ~ "exercise_no_controls",
      TRUE ~ "exercise_with_controls"
    ),
    contrast_type = factor(contrast_type,
                           levels = c("exercise_with_controls",
                                      "exercise_no_controls",
                                      "Endur_vs_Resist",
                                      "baseline",
                                      "control_only")),
    contrast_order = 1:n(),
    contrast_category = sub("([^.]+)\\.[^-]+ - ([^.]+).*",
                            "\\1-\\2",
                            contrast_short),
    contrast_category = gsub("Control", "CON", contrast_category),
    contrast_category = gsub("Endur", "EE", contrast_category),
    contrast_category = gsub("Resist", "RE", contrast_category),
    across(.cols = c(contrast, contrast_short,
                     contrast_type, contrast_category),
           .fns = ~ factor(.x, levels = unique(.x)))
  ) %>%
  # contrast_left is the left-hand group.timepoint of the contrast, e.g. "Endur.post_10_min".
  # str_split_fixed(..., 2) splits on the FIRST separator only, which is load-bearing for the
  # timepoint: "post_3.5_4_hr" itself contains a "." and must not be split further.
  dplyr::mutate(contrast_left = stringr::str_split_fixed(contrast_short, " - ", 2)[,1],
                randomGroupCode = case_when(
                  contrast_category == "Endur_vs_Resist" ~ "ADUEndur - ADUResist",
                  TRUE ~ paste0("ADU", stringr::str_split_fixed(contrast_left, "\\.", 2)[,1])
                ),
                Timepoint = as.factor(stringr::str_split_fixed(contrast_left, "\\.", 2)[,2])) %>%
  dplyr::select(-contrast_left) %>%
  relocate(contrast_order, .before = everything()) %>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("pre_exercise",
                                       "during_20_min",
                                       "during_40_min",
                                       "post_10_min",
                                       "post_15_30_45_min",
                                       "post_3.5_4_hr",
                                       "post_24_hr")))
  #make sure the timepoints are in order of actual timepoints, instead of otherwise

setDT(CONTRAST_CONVERTER)
save(CONTRAST_CONVERTER, file = file.path(out_dir, "CONTRAST_CONVERTER.rda"),
     compress = TRUE, version = 3)
message(sprintf("CONTRAST_CONVERTER: %d contrast(s), %d column(s) -> %s",
                nrow(CONTRAST_CONVERTER), ncol(CONTRAST_CONVERTER), out_dir))

CONTRAST_LEVELS <- levels(CONTRAST_CONVERTER[["contrast"]])

# ---- the DA freeze inventory ------------------------------------------------------------
# Anchored on `_da_<model>_v<x.y>.txt` so both model families match (dream for everything
# fit here, malax-glmm for the vendored methylcap tier) while the metadata_*, qc-norm_* and
# imputed_* files beside them do not. Same pattern step 11 settled on.
freeze_dir <- file.path(ROOT, "staging", "freeze")
da_files <- list.files(freeze_dir, pattern = "_da_.*_v[0-9.]+\\.txt$",
                       recursive = TRUE, full.names = TRUE)
if (length(da_files) == 0)
  stop("no DA freeze files under staging/freeze/*/da/ — run step 09 first")

OME_TISSUE_CODE <- readRDS(file.path(ROOT, "scripts", "00_preflight", "data", "OME_TISSUE_CODE.rds"))

# METABOLOMICS_CVS (step 05) — the CV-based metabolite dedup table. See collapse_metab_cv().
local({
  f <- file.path(out_dir, "METABOLOMICS_CVS.rda")
  if (!file.exists(f))
    stop("missing METABOLOMICS_CVS.rda (step 05) — needed to collapse the metab DA feature space: ", f)
  e <- new.env(); load(f, envir = e); METABOLOMICS_CVS <<- e$METABOLOMICS_CVS
})

# human-precovid-sed-adu_<tissue_code>_<ome>_da_<model>-<type>_v<x.y>.txt
parse_da_name <- function(f) {
  tok <- sub("_da_.*$", "", sub("^human-precovid-sed-adu_", "", basename(f)))
  list(tissue_code = sub("_.*$", "", tok), ome = sub("^[^_]+_", "", tok))
}

inv <- lapply(da_files, function(f) {
  p <- parse_da_name(f)
  data.frame(file = f, tissue_code = p$tissue_code, ome = p$ome, stringsAsFactors = FALSE)
})
inv <- do.call(rbind, inv)

# tissue from the BIC code
inv$tissue <- vapply(seq_len(nrow(inv)), function(i) {
  hit <- OME_TISSUE_CODE[OME_TISSUE_CODE$ome == inv$ome[i] &
                           OME_TISSUE_CODE$tissue_code == inv$tissue_code[i], ]
  if (nrow(hit) == 0)
    stop("cannot map ", inv$tissue_code[i], "/", inv$ome[i],
         " back to a tissue via OME_TISSUE_CODE: ", basename(inv$file[i]))
  as.character(hit[["tissue"]][1])
}, character(1))

# The object key. Every metabolomics platform except the two clinical ones collapses to a
# single "metab" object per tissue, carrying the platform in its own column — upstream's
# bind_rows(.id = "platform"). metab-t-clinical stays on its own key (see the header), as do
# prot-* and transcript-rna-seq, which are one object per fitted table.
inv$ome_group <- ifelse(grepl("^metab-", inv$ome) & inv$ome != "metab-t-clinical",
                        "metab", inv$ome)

# {TISSUE}_{OME}_DA. transcript-rna-seq is the one ome whose object name is not its assay
# code spelled out; upstream abbreviates it to TRNSCRPT and downstream code uses that name.
ome_abbr <- function(ome) if (ome == "transcript-rna-seq") "TRNSCRPT" else toupper(gsub("-", "_", ome))
inv$object <- sprintf("%s_%s_DA", toupper(inv$tissue), vapply(inv$ome_group, ome_abbr, character(1)))

message(sprintf("DA freeze: %d file(s) -> %d object(s)",
                nrow(inv), length(unique(inv$object))))

# ---- read + reshape ---------------------------------------------------------------------
read_da <- function(f) {
  # fread, not read.csv: SEVERAL VERY LARGE MATRICES AT ONCE. Unlike step 11 this step keeps
  # every column, so no subset is possible — the non-epigen DA tier is ~890 MB across 43
  # files, the blood transcriptomics table alone 281 MB. read.csv on that set is minutes of
  # wall clock and several GB of peak memory; fread reads it in one pass.
  data.table::fread(f, sep = "\t", header = TRUE, data.table = TRUE)
}

# Collapse one metab platform's DA table onto the RefMet/lowest-CV feature space.
#
# THIS IS THE DEDUP STEP 09 DEFERRED. generate_metab_da.R fits with
# `remove_redundant_metab = FALSE`, overriding load_qc_local()'s default, and says why:
# "the CV-based dedup belongs downstream of the DA, not in front of it". Downstream is here.
# So the metab freeze is deliberately keyed on NATIVE metabolite names over the full
# un-deduped feature space ("CoA(4:0)", "Succinate", "DG(34:3)"), and this step is what
# reduces it to the released feature space, keyed on RefMet names ("Crotonoyl-CoA",
# "Succinic acid").
#
# Two things happen, in the same order and with the same semantics as load_qc_local():
#   1. keep only features flagged lowest_CV == "yes" for this tissue x platform — the same
#      analyte is measured at several sites, and only the least-variable site's copy is kept
#   2. rename feature_id to refmet_name where one exists (4.4% have none; those keep the
#      native id, which is load_qc_local()'s case_when fallback)
#
# Getting this wrong is not cosmetic. Without it the DA and the *_QC objects sit in two
# different namespaces, and step 11 — which restricts each qc-norm matrix to the features
# present in its DA — silently intersects to nothing: 9 of 47 tissue x ome combos came out
# empty, and pivot_longer() then failed on a zero-column frame rather than reporting it.
#
# metab-t-clinical is NOT collapsed, and never reaches this function: load_qc_local() excludes
# it from the same dedup (`platforms != "metab-t-clinical"`), and it is its own ome_group here
# rather than part of the metab stack, so the two agree by construction.
#
# The lookup is unique per (tissue, assay, feature_id) among lowest_CV == "yes" rows, so the
# match cannot multiply DA rows; that is asserted rather than assumed.
collapse_metab_cv <- function(dt, tissue, platform) {
  lk <- METABOLOMICS_CVS[METABOLOMICS_CVS$tissue == tissue &
                           METABOLOMICS_CVS$assay == platform &
                           METABOLOMICS_CVS$lowest_CV == "yes",
                         c("feature_id", "refmet_name")]
  lk <- unique(lk)
  before <- length(unique(dt$feature_id))

  # A platform can lose every feature to the dedup: each of its analytes is measured at a
  # lower CV somewhere else, so nothing it carries is flagged lowest_CV == "yes". That is a
  # legitimate result of the filter and not a malformed lookup, so it is reported and the
  # platform contributes no rows to the stack. blood/metab-t-imm-crt is the case that occurs
  # (see the header). The zero-row table keeps the freeze's column types so the platform still
  # stacks cleanly in rbindlist().
  if (nrow(lk) == 0) {
    message(sprintf("    CV/RefMet collapse: %s/%s has no lowest_CV == \"yes\" feature(s) — all %d feature(s) drop, platform contributes no rows",
                    tissue, platform, before))
    out <- dt[0L]
    attr(out, "collapse_note") <- sprintf("%s: %d -> 0 feature(s) [none lowest_CV]",
                                          platform, before)
    return(out)
  }
  if (anyDuplicated(lk$feature_id))
    stop(sprintf("%s/%s: METABOLOMICS_CVS maps one feature_id to >1 refmet_name", tissue, platform))

  # `lowest_CV == "yes"`, not `!= "no"`: a feature absent from METABOLOMICS_CVS joins to NA,
  # and load_qc_local()'s filter(!lowest_CV == "no") drops NA too (dplyr::filter drops NA).
  # Both sides therefore keep exactly the features the table flags "yes".
  dt <- dt[dt$feature_id %in% lk$feature_id]
  # as.character on refmet_name is load-bearing: METABOLOMICS_CVS stores it as a FACTOR, and
  # ifelse() drops attributes from its branches — a factor arm comes back as its integer
  # codes. Left as a factor this silently renamed every metabolite to a number ("79", "182"),
  # which still passed a "no native ids survived" check and only surfaced downstream as step
  # 11 intersecting to zero features.
  map <- stats::setNames(as.character(lk$refmet_name), as.character(lk$feature_id))
  new <- unname(map[dt$feature_id])
  dt[, feature_id := ifelse(!is.na(new) & nzchar(new), new, feature_id)]
  attr(dt, "collapse_note") <- sprintf("%s: %d -> %d feature(s)", platform, before,
                                       length(unique(dt$feature_id)))
  dt
}

# Port of MotrpacHumanPreSuspensionAnalysis:::.process_raw_DA() for a single table. Upstream
# takes the whole named list and derives tissue/assay by splitting each "<tissue>.<ome>"
# name; the split happens in the caller here, so they arrive as arguments instead.
process_da <- function(dt, tissue, assay) {
  # contrast is the join key back to the converter. It is character in the freeze and a
  # factor in CONTRAST_CONVERTER, so it is cast for the merge and re-levelled after —
  # matching on the factor's integer codes would silently mis-key.
  dt[, contrast := as.character(contrast)]
  cc <- data.table::copy(CONTRAST_CONVERTER)
  cc[, contrast := as.character(contrast)]

  unknown <- setdiff(unique(dt$contrast), cc$contrast)
  if (length(unknown))
    stop(sprintf("%s/%s: %d contrast(s) absent from contrast_converter.txt, e.g. %s",
                 tissue, assay, length(unknown), paste(utils::head(unknown, 2), collapse = " | ")))

  dt <- merge(x = dt, y = cc, by = "contrast", all.x = TRUE, all.y = FALSE)

  dt[, `:=`(tissue = tissue, assay = assay)]
  dt[, `:=`(
    contrast   = factor(x = contrast, levels = CONTRAST_LEVELS),
    feature_id = as.character(feature_id),
    full_model = as.factor(full_model)
  )]
  dt[, contrast_order := NULL]
  if ("platform" %in% colnames(dt)) dt[, platform := as.factor(platform)]
  dt[, contrast := droplevels(contrast)]

  # Column order: identity block, then the converter's columns, then whatever statistics the
  # table happens to carry, in the order the freeze wrote them. intersect() throughout, so a
  # table missing a column (methylcap has no z.std) reshapes rather than errors.
  new_order <- intersect(
    x = c("tissue", "assay", "platform", "full_model", colnames(CONTRAST_CONVERTER)),
    y = colnames(dt)
  )
  data.table::setcolorder(x = dt, neworder = new_order)
  if ("z.std" %in% colnames(dt))
    data.table::setcolorder(x = dt, neworder = "z.std", before = "p_value")

  keys <- intersect(x = c("full_model", "contrast", "platform", "p_value", "feature_id"),
                    y = colnames(dt))
  data.table::setkeyv(x = dt, cols = keys)
  dt[]
}

written <- character(0)
for (obj in sort(unique(inv$object))) {
  rows <- inv[inv$object == obj, , drop = FALSE]
  tissue <- unique(rows$tissue)
  if (length(tissue) != 1) stop(obj, ": spans >1 tissue: ", paste(tissue, collapse = ", "))

  parts <- lapply(seq_len(nrow(rows)), function(i) read_da(rows$file[i]))
  names(parts) <- rows$ome

  if (unique(rows$ome_group) == "metab") {
    # Upstream stacks the platforms with bind_rows(.id = "platform"), taking the platform
    # from the list name. The freeze already writes a platform column (assay = "metab",
    # platform = the ome), so the name is verified against it instead of overwriting it —
    # if the two ever disagree the file is mislabelled and that should surface here.
    for (i in seq_along(parts)) {
      if ("platform" %in% colnames(parts[[i]])) {
        got <- unique(as.character(parts[[i]]$platform))
        if (length(got) != 1 || got != names(parts)[i])
          stop(sprintf("%s: platform column [%s] disagrees with filename ome [%s]",
                       basename(rows$file[i]), paste(got, collapse = ", "), names(parts)[i]))
      } else {
        parts[[i]][, platform := names(parts)[i]]
      }
    }
    # Collapse each platform onto the RefMet/lowest-CV space before stacking. Per platform,
    # because the METABOLOMICS_CVS lookup is keyed on (tissue, assay = platform, feature_id).
    notes <- character(0)
    for (i in seq_along(parts)) {
      parts[[i]] <- collapse_metab_cv(parts[[i]], tissue = tissue, platform = names(parts)[i])
      notes <- c(notes, attr(parts[[i]], "collapse_note"))
    }
    message(sprintf("    CV/RefMet collapse: %s", paste(notes, collapse = "; ")))
  }

  dt <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  # assay is the object's ome, not whatever the freeze wrote: for the metab stack the freeze
  # says "metab" per file, and .process_raw_DA() likewise overwrites assay from the list name.
  dt <- process_da(dt, tissue = tissue, assay = unique(rows$ome_group))

  assign(obj, dt)
  save(list = obj, file = file.path(out_dir, paste0(obj, ".rda")), compress = TRUE, version = 3)
  written <- c(written, obj)
  message(sprintf("  %-28s %2d file(s) %8d rows x %2d cols  %s",
                  obj, nrow(rows), nrow(dt), ncol(dt),
                  if ("platform" %in% colnames(dt))
                    sprintf("%d platform(s)", length(unique(as.character(dt$platform)))) else ""))
  rm(list = obj); rm(dt, parts); invisible(gc(verbose = FALSE))
}

message(sprintf("DA_ASSEMBLE: wrote %d object(s) + CONTRAST_CONVERTER -> %s",
                length(written), out_dir))
