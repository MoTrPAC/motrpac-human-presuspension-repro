#!/usr/bin/env Rscript
# Leaf data object: OME_TISSUE_CODE  (deps=-, non-gated)
# Adapted from .../data-raw/ome_tissue_object.R. Instead of a gated `gsutil ls -R`
# on the results bucket, this derives the (tissue, tissue_code, ome) map from the
# QUANTID_BUCKET_FILES data object (../data/QUANTID_BUCKET_FILES.rds), so it needs
# no consortium access. `tissue` is the overarching category (blood/muscle/adipose)
# from the vendored .find_tissue(); `ome` is the assay_class. All assays are kept,
# INCLUDING clinical-chemistry (lab-*); only qa-qc rows (null tissue_code) drop.
suppressWarnings(suppressMessages(library(dplyr)))
source(file.path(Sys.getenv("PREFLIGHT_DR"), "lib", "helpers.R"))

rds <- file.path(.out_dir(), "QUANTID_BUCKET_FILES.rds")
if (!file.exists(rds)) stop("missing QUANTID_BUCKET_FILES.rds — build it first: ", rds)
files <- readRDS(rds)

keep <- (files$is_qa_qc %in% FALSE) & !is.na(files$tissue_code)
# exclude the t10 ATAC reference standard (t10-muscle-powder epigen-atac-seq) —
# matches upstream ome_tissue_object.R, which drops "t10 atac" as reference only.
keep <- keep & !(grepl("t10", files$tissue_code) & files$assay_class == "epigen-atac-seq")
sub <- files[keep, c("tissue_code", "assay_class")]

OME_TISSUE_CODE <- sub %>%
  transmute(tissue = vapply(tissue_code, function(tc) {
              t <- .find_tissue(tc); if (is.null(t)) NA_character_ else t
            }, character(1)),
            tissue_code = tissue_code,
            ome = assay_class) %>%
  filter(!is.na(tissue)) %>%
  distinct(.keep_all = TRUE)
OME_TISSUE_CODE <- as.data.frame(OME_TISSUE_CODE)

# manually add methylseq + clinical rows absent from the quant-id tier, row by
# row. distinct() dedupes any overlap.
OME_TISSUE_CODE <- rbind(OME_TISSUE_CODE, c("blood", "t03-edta", "epigen-methylcap-seq"))
OME_TISSUE_CODE <- rbind(OME_TISSUE_CODE, c("muscle", "t06-muscle", "epigen-methylcap-seq"))
OME_TISSUE_CODE <- rbind(OME_TISSUE_CODE, c("adipose", "t11-adipose", "epigen-methylcap-seq"))
# clinical chemistry is curated upstream and written as these omes
OME_TISSUE_CODE <- rbind(OME_TISSUE_CODE, c("blood", "t02-plasma", "metab-t-clinical"))
OME_TISSUE_CODE <- rbind(OME_TISSUE_CODE, c("blood", "t02-plasma", "prot-clinical"))
OME_TISSUE_CODE <- dplyr::distinct(OME_TISSUE_CODE, .keep_all = TRUE)

saveRDS(OME_TISSUE_CODE, file.path(.out_dir(), "OME_TISSUE_CODE.rds"))
message("OME_TISSUE_CODE: ", nrow(OME_TISSUE_CODE), " rows (from quant-id catalog, incl clinical-chemistry) -> ", .out_dir())
