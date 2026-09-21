#!/usr/bin/env Rscript
# Stage 1 data objects: the *_QC objects  (deps: QC_NORM, pheno)
# Adapted from MotrpacHumanPreSuspensionData/data-raw/qc_norm_results.R.
#
# The upstream script calls the package-internal .load_qc(), which lists a staging
# bucket and dl_read_gcp's every freeze file into a nested [[tissue]][[ome]] list
# (each leaf a list of qc_norm / feature_metadata / sample_metadata), then renames
# to {TISSUE}_{OME}_QC and use_data()'s each. This adaptation instead reads the
# LOCAL metadata + qc-norm files that the generate_normalized_expression stems
# wrote into staging/freeze/*/, and save()s each object as a local .rda in
# scripts/10_build_data/data/.
#
# Two deliberate deviations from the upstream script:
#   - Imputed prot platforms are NOT built: `imputed` freeze files are skipped, so
#     the prot-pr/ph *_QC objects carry no qc_imputed component.
#   - Clinical omes ARE included: prot-clinical and metab-t-clinical are packaged
#     as BLOOD_PROT_CLINICAL_QC / BLOOD_METAB_T_CLINICAL_QC (the upstream ome
#     factor omitted them).
#
# Blood transcriptomics is split into blood_transcript_1 / blood_transcript_2
# (qc_norm columns 1:450 / 451:end) to stay under GitHub's file-size cap, matching
# how load_qc() cbind-recombines them.
#
# NOTE: these data/*.rda are staging — they will eventually be organized as the
# canonical *_QC data objects in the MotrpacHumanPreSuspensionData package.
suppressWarnings(suppressMessages({ library(dplyr); library(data.table) }))

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
# qc_helpers.R provides load_pheno() and loads the `pheno` + OME_TISSUE_CODE globals.
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

out_dir    <- file.path(ROOT, "scripts", "10_build_data", "data")
freeze_dir <- file.path(ROOT, "staging", "freeze")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pheno_all <- load_pheno(load_acute_only = FALSE)  # upstream uses load_acute_only = FALSE

# ome x tissue combos that become *_QC objects: everything catalogued in
# OME_TISSUE_CODE except the clinical chemistry labs. Clinical
# proteomics/metabolomics ARE kept, and so is epigenetics — see the methylcap
# note below.
combos <- OME_TISSUE_CODE %>%
  dplyr::filter(!grepl("^lab-", ome)) %>%
  dplyr::mutate(token = paste0(tissue_code, "_", ome))

ome_label <- function(ome) ifelse(ome == "transcript-rna-seq", "TRNSCRPT",
                                  toupper(gsub("-", "_", ome)))

# --- assemble each leaf from the local freeze files -------------------------
# read.csv (not fread): fread mis-detects the last row of some large qc-norm
# matrices as a footer (files have no trailing newline), silently dropping a
# feature. quote/comment disabled so BIC fields are taken literally, and
# check.names = FALSE keeps the vial-label column names intact.
read_freeze <- function(f) utils::read.csv(f, sep = "\t", header = TRUE, check.names = FALSE,
                                             stringsAsFactors = FALSE, quote = "", comment.char = "")

# Each (token, type) must resolve to exactly one freeze file. Multiple matches mean a
# stale prior-version file is lingering (e.g. _v1.2 beside _v1.4) — fail loudly so it
# gets cleaned up rather than silently picking one. 0 matches -> NULL (skip).
.one_file <- function(files, what, token) {
  if (length(files) == 0) return(NULL)
  if (length(files) == 1) return(files)
  stop(sprintf("multiple %s freeze files for %s (remove stale versions):\n  %s",
               what, token, paste(basename(files), collapse = "\n  ")))
}

freeze_files <- list.files(freeze_dir, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
freeze_files <- freeze_files[!grepl("imputed", freeze_files)]   # skip imputed prot platforms

qc_list <- list()
for (i in seq_len(nrow(combos))) {
  tk <- combos$token[i]
  files_i <- freeze_files[grepl(tk, basename(freeze_files), fixed = TRUE)]
  if (length(files_i) == 0) next
  obj_name <- paste0(toupper(combos$tissue[i]), "_", ome_label(combos$ome[i]), "_QC")
  leaf <- list()

  qn <- .one_file(files_i[grepl("qc-norm", files_i)], "qc-norm", tk)
  if (!is.null(qn)) {
    m <- read_freeze(qn)
    rownames(m) <- m[, 1]; m[, 1] <- NULL
    leaf$qc_norm <- m
  }
  fm <- .one_file(files_i[grepl("metadata_features", files_i)], "metadata_features", tk)
  if (!is.null(fm)) leaf$feature_metadata <- read_freeze(fm)

  sm <- .one_file(files_i[grepl("metadata_samples", files_i)], "metadata_samples", tk)
  if (!is.null(sm)) {
    s <- read_freeze(sm)
    # merge with pheno on vialLabel, dropping columns shared with pheno_data
    # (keep vialLabel; drop=FALSE guards clinical omes whose sample metadata
    # shares every column with pheno_data).
    shared <- setdiff(intersect(names(s), names(pheno_all$pheno_data)), "vialLabel")
    s <- s[, !(names(s) %in% shared), drop = FALSE]
    leaf$sample_metadata <- merge(s, pheno_all$pheno_data, by = "vialLabel")
  }
  # keep the matrix and its metadata in step: pheno covers only the sed-adult cohort, so a
  # column with no pheno row would otherwise survive with no covariates. This matters for
  # methylcap, whose beta-value matrix is copied into the freeze un-subset.
  if (!is.null(leaf$qc_norm) && !is.null(leaf$sample_metadata))
    leaf$qc_norm <- leaf$qc_norm[, colnames(leaf$qc_norm) %in%
                                   as.character(leaf$sample_metadata$vialLabel), drop = FALSE]
  qc_list[[obj_name]] <- leaf
}

# --- blood transcriptomics split --------------------------------------------
if (!is.null(qc_list[["BLOOD_TRNSCRPT_QC"]])) {
  bt <- qc_list[["BLOOD_TRNSCRPT_QC"]]
  ncol_bt <- ncol(bt$qc_norm)
  if (is.null(ncol_bt) || ncol_bt <= 450)
    stop("BLOOD_TRNSCRPT qc_norm has ", ncol_bt, " columns; expected > 450 to split")
  blood_transcript_1 <- bt; blood_transcript_1$qc_norm <- bt$qc_norm[, 1:450]
  blood_transcript_2 <- bt; blood_transcript_2$qc_norm <- bt$qc_norm[, 451:ncol_bt]
  qc_list[["BLOOD_TRNSCRPT_QC"]] <- NULL
  qc_list[["blood_transcript_1"]] <- blood_transcript_1
  qc_list[["blood_transcript_2"]] <- blood_transcript_2
}

# --- save each object as a local .rda ---------------------------------------
for (nm in names(qc_list)) {
  assign(nm, qc_list[[nm]])
  save(list = nm, file = file.path(out_dir, paste0(nm, ".rda")), compress = TRUE, version = 3)
}
message("qc_norm_results: wrote ", length(qc_list), " objects -> ", out_dir)
message("  ", paste(names(qc_list), collapse = ", "))
