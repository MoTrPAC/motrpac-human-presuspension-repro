#!/usr/bin/env Rscript
# Step 09 DA stem: clinical chemistry (blood; prot-clinical and metab-t-clinical).
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_differential_analysis/
# generate_clinical_DA.R.
#
# Output: staging/freeze/{proteomics,metabolomics-targeted}/da/
#   human-precovid-sed-adu_t02-plasma_{prot-clinical,metab-t-clinical}_da_dream-acute_v<x.y>.txt

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "09_build_da", "da_common.R"))

# ---- Upstream input assembly, verbatim from .../generate_clinical_DA.R ----
#
# clin_chemistry = MotrpacHumanPreSuspensionData::load_clinical_data()[["chemistry"]]
# clin_chem_df = clin_chemistry %>%
#   dplyr::bind_rows() %>%
#   tibble::column_to_rownames("analyte_name")
# #just plug in a random targeted metab platform for metadata and covariates
# metab_metadata = MotrpacHumanPreSuspensionData::load_pheno()[["pheno_data"]]
# shared_samples = intersect(metab_metadata$vialLabel, colnames(clin_chem_df))
# metab_metadata = metab_metadata %>% dplyr::filter(vialLabel %in% shared_samples)
# clin_chem_df = clin_chem_df %>%
#   dplyr::select(dplyr::any_of(shared_samples)) %>%
#   dplyr::mutate(dplyr::across(dplyr::everything(), ~log2(.x)))
#


CLINICAL_OMES <- c("prot-clinical", "metab-t-clinical")

# PARALLEL PATH: `parallel` reaches run_dream()'s MulticoreParam fork backend. See
# VARIANCEPARTITION_PARALLEL_CORES in config/pipeline.env for why it is fork and not SOCK. This
# stem is the sharpest case: every prot-clinical analyte carries NAs, so under the parallel
# backend every model fails rather than a few going missing.
fit_clinical_ome = function(desired_ome, tissue = "blood", parallel = da_parallel_enabled()){
  qc = load_qc_local(selected_tissues = tissue,
                     selected_omes = desired_ome,
                     load_acute_only = FALSE,
                     load_clinical = TRUE)
  leaf = qc[[tissue]][[desired_ome]]
  if (is.null(leaf))
    stop("no local *_QC object for ", tissue, "/", desired_ome, " — run step 08 first")

  metadata = leaf[["sample_metadata"]] %>% dplyr::filter(visitcode == "ADU_BAS")
  rownames(metadata) = metadata$vialLabel

  # absolute values in, log2 out — the clinical qc-norm stem deliberately does not transform
  data_matrix = leaf[["qc_norm"]] %>%
    dplyr::select(as.character(metadata$vialLabel)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ log2(.x)))

  process_metadata = process_covariates(meta = metadata,
                                        selected_ome = desired_ome,
                                        tissue_input = tissue,
                                        include_technical = TRUE)

  # run_dream() rather than run_da_models(): the latter writes the freeze file immediately,
  # and the BH correction below spans both omes, so writing has to wait until both are fit.
  # Sample alignment is explicit inside run_dream(), which reorders the covariate rows with
  # match(colnames(expression_object), rownames(meta_matrix)) before fitting.
  fit = run_dream(expression_object = data_matrix,
                  model_type = DA_MODEL_TYPE,
                  process_metadata = process_metadata,
                  voom = FALSE,
                  parallel = parallel)

  out = .convert_dream_output(fit,
                              metadata = process_metadata$original_meta,
                              tissue = tissue,
                              formula = process_metadata[["full_formula"]],
                              ome = desired_ome)
  message(sprintf("%s / %s: %d analyte(s), %d rows", tissue, desired_ome,
                  nrow(data_matrix), nrow(out)))
  return(out)
}

# ---- multiple-testing correction ----
# BH adjustment across all 9 features.
clin_all <- dplyr::bind_rows(lapply(CLINICAL_OMES, fit_clinical_ome))
clin_all$adj_p_value <- stats::ave(clin_all$p_value, clin_all$contrast,
                                   FUN = function(p) stats::p.adjust(p, method = "BH"))
message(sprintf("BH within contrast across %d analyte(s) in %d contrast(s)",
                dplyr::n_distinct(clin_all$feature_id), dplyr::n_distinct(clin_all$contrast)))

rows_written <- 0
for (desired_ome in CLINICAL_OMES) {
  ome_da = clin_all %>% dplyr::filter(assay == desired_ome)

  # .convert_dream_output() sets assay = ome for every stem, and the filter above splits the
  # two omes on that value. The metab side then takes the identity the other metab platforms
  # use (generate_metab_da.R): assay = "metab" — the value HUMAN_FEATURE_TO_GENE keys on —
  # with the ome in its own platform column. The prot platforms name themselves in assay, so
  # prot-clinical is written as fit.
  if (grepl("^metab", desired_ome))
    ome_da = ome_da %>%
      dplyr::mutate(platform = assay, assay = "metab") %>%
      dplyr::relocate(platform, .after = assay)

  da_path = file.path(.STAGING, "freeze", .freeze_subdir_for_ome(desired_ome), "da")
  dir.create(da_path, recursive = TRUE, showWarnings = FALSE)
  write_with_path_name(ome_da,
                       local_path = da_path,
                       ome = desired_ome,
                       tissue = "blood",
                       data_category = "da",
                       data_details = paste0("dream-", DA_MODEL_TYPE))
  message(sprintf("blood / %s: %d rows across %d contrast(s) -> %s",
                  desired_ome, nrow(ome_da), dplyr::n_distinct(ome_da$contrast), da_path))
  rows_written <- rows_written + nrow(ome_da)
}

# the split must be lossless -- every row lands in exactly one ome
stopifnot(rows_written == nrow(clin_all))
