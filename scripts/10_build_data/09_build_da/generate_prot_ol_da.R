#!/usr/bin/env Rscript
# Step 09 DA stem: prot-ol (Olink targeted proteomics, blood only).
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_differential_analysis/
# generate_prot_ol_DA.R.
#
# Technical covariates are excluded here (include_technical = FALSE) because the Olink
# assay's own normalization already accounts for the major technical variation, and no
# voom weighting is applied — the model is fit directly on normalized expression.
#
# Output: staging/freeze/proteomics/da/
#   human-precovid-sed-adu_t02-plasma_prot-ol_da_dream-acute_v<x.y>.txt

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "09_build_da", "da_common.R"))

# ---- Upstream build, verbatim from .../generate_differential_analysis/generate_prot_ol_DA.R ----
#
# .generate_prot_ol_inputs = function(repo_local_dir,
#                                     model_type,
#                                     tissue,
#                                     parallel = F){
#   desired_ome = 'prot-ol'
#
#   da_path = file.path(repo_local_dir, "data", "tmp", "freeze_DA/") #path for output
#   dir.create(da_path, recursive = TRUE, showWarnings = FALSE)
#
#   prot_ol_data = load_qc(selected_tissues=tissue,
#                          selected_omes=desired_ome,
#                          load_acute_only=FALSE)
#   if (length(prot_ol_data) > 0) {
#     metadata = prot_ol_data[[tissue]][[desired_ome]][['sample_metadata']]
#     if (model_type == "acute") metadata = metadata %>% dplyr::filter(visitcode == 'ADU_BAS')  #filter just to the initial acute bout
#     rownames(metadata) = metadata$vialLabel
#
#     data_matrix = prot_ol_data[[tissue]][[desired_ome]][['qc_norm']] %>%
#       dplyr::select(as.character(metadata$vialLabel))
#
#     process_metadata = process_covariates(meta = metadata,
#                                           selected_ome = desired_ome,
#                                           tissue_input = tissue,
#                                           include_technical = F)
#
#     .run_models(repo_local_dir = repo_local_dir,
#                 model_type = model_type,
#                 expression_object = data_matrix,
#                 process_metadata = process_metadata,
#                 tissue = tissue,
#                 ome = desired_ome,
#                 voom = FALSE,
#                 parallel = parallel)
#   }
# }
#
# Adapted below: load_qc() (installed MotrpacHumanPreSuspensionData) becomes load_qc_local()
# reading this stage's step-08 *_QC.rda, and .run_models() becomes run_da_models() writing
# into staging/freeze. The modeling itself is unchanged.

# PARALLEL PATH: `parallel` reaches run_dream()'s MulticoreParam fork backend. See
# VARIANCEPARTITION_PARALLEL_CORES in config/pipeline.env for why it is fork and not SOCK.
generate_prot_ol_da = function(tissue = "blood", parallel = da_parallel_enabled()){
  desired_ome = "prot-ol"

  prot_ol_data = load_qc_local(selected_tissues = tissue,
                               selected_omes = desired_ome,
                               load_acute_only = FALSE)
  if (length(prot_ol_data) == 0)
    stop("no local *_QC object for ", tissue, "/", desired_ome, " — run step 08 first")

  metadata = prot_ol_data[[tissue]][[desired_ome]][["sample_metadata"]]
  metadata = metadata %>% dplyr::filter(visitcode == "ADU_BAS")   # initial acute bout only
  rownames(metadata) = metadata$vialLabel

  data_matrix = prot_ol_data[[tissue]][[desired_ome]][["qc_norm"]] %>%
    dplyr::select(as.character(metadata$vialLabel))

  process_metadata = process_covariates(meta = metadata,
                                        selected_ome = desired_ome,
                                        tissue_input = tissue,
                                        include_technical = FALSE)

  # Sample alignment is explicit inside run_da_models(): run_dream() reorders the covariate
  # rows with match(colnames(expression_object), rownames(meta_matrix)) before fitting, so
  # the design always follows the expression columns rather than their incoming order.
  out = run_da_models(expression_object = data_matrix,
                      process_metadata = process_metadata,
                      tissue = tissue,
                      ome = desired_ome,
                      voom = FALSE,
                      parallel = parallel)
  return(invisible(out))
}

generate_prot_ol_da()
