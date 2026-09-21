#!/usr/bin/env Rscript
# Leaf data objects: figure color/abbreviation palettes  (deps=-, non-gated)
# Vendored verbatim from .../data-raw/COLORS_ABBREVIATIONS.R, with
# usethis::use_data(X) replaced by saveRDS(X, <data>/X.rds) so it builds outside
# an R package. Fully self-contained (hardcoded literals; no external input).
source(file.path(Sys.getenv("PREFLIGHT_DR"), "lib", "helpers.R"))
.OUT <- .out_dir()
.save <- function(obj, name) saveRDS(obj, file.path(.OUT, paste0(name, ".rds")))

HUMAN_SEX_COLORS = c(male = "#5555ff", Male = "#5555ff",
                     female = "#f95c6f", Female = "#f95c6f")
.save(HUMAN_SEX_COLORS, "HUMAN_SEX_COLORS")

HUMAN_TISSUE_COLORS = c(blood = "#d7191c", muscle = "#abd9e9", adipose = "#ffffbf")
.save(HUMAN_TISSUE_COLORS, "HUMAN_TISSUE_COLORS")

HUMAN_TISSUE_ABBR = c(blood = "#BLO", muscle = "#MUS", adipose = "#ADI")
.save(HUMAN_TISSUE_ABBR, "HUMAN_TISSUE_ABBR")

HUMAN_EXERCISE_GROUP_COLORS = c(ADUResist = "#1b9e77", ADUEndur = "#d95f02",
                                ADUControl = "#7570b3")
.save(HUMAN_EXERCISE_GROUP_COLORS, "HUMAN_EXERCISE_GROUP_COLORS")

HUMAN_ACUTE_TIMEPOINT_COLORS = c(pre_exercise = "#BEBEBE", Pre_exercise = "#BEBEBE",
                                 During = "#FDE725", during_20_min = "#FDE725",
                                 during_40_min = "#BAD071", post_10_min = "#D1BBD7",
                                 Early = "#AE76A3", post_15_30_45_min = "#AE76A3",
                                 post_15_min = "#AE76A3", post_30_min = "#AE76A3",
                                 post_45_min = "#AE76A3", Mid = "#882E72",
                                 post_3.5_4_hr = "#882E72", post_4_hr = "#882E72",
                                 Late = "#61194F", post_24_hr = "#61194F")
.save(HUMAN_ACUTE_TIMEPOINT_COLORS, "HUMAN_ACUTE_TIMEPOINT_COLORS")

HUMAN_OME_COLORS = c(Transcriptomics = "#377EB8", "transcript-rna-seq" = "#377EB8",
                     Proteomics = "#228833", "Proteomics (MS)" = "#228833",
                     "prot-pr" = "#228833", "Proteomics (Olink)" = "#228833",
                     "prot-ol" = "#228833", Metabolomics = "#6D4B08",
                     "metab-u-hilicpos" = "#6D4B08", "metab-u-ionpneg" = "#6D4B08",
                     "metab-u-lrpneg" = "#6D4B08", "metab-u-lrppos" = "#6D4B08",
                     "metab-u-rpneg" = "#6D4B08", "metab-u-rppos" = "#6D4B08",
                     "metab-t-amines" = "#6D4B08", "metab-t-conv" = "#6D4B08",
                     "metab-t-imm-crt" = "#6D4B08", "metab-t-imm-glc" = "#6D4B08",
                     "metab-t-imm-ins" = "#6D4B08", "metab-t-oxylipneg" = "#6D4B08",
                     "metab-t-tca" = "#6D4B08", "metab-t-nuc" = "#6D4B08",
                     "metab-t-acoa" = "#6D4B08", "metab-t-ka" = "#6D4B08",
                     "epigen-atac-seq" = "#882255",
                     ATAC = "#882255", `Chromatin Accessibility (ATAC)` = "#882255",
                     Phosphoproteomics = "#F3A02B", "prot-ph" = "#F3A02B",
                     "epigen-methylcap-seq" = "#D687B5", Methylation = "#D687B5")
.save(HUMAN_OME_COLORS, "HUMAN_OME_COLORS")

message("COLORS_ABBREVIATIONS: 6 objects -> ", .OUT)
