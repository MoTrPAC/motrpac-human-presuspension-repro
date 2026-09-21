#!/usr/bin/env Rscript
# Step 09 DA stem: prot-ph (phosphoproteomics; adipose + muscle).
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_differential_analysis/
# generate_prot_ph_DA.R.
#
# Structurally identical to the prot-pr stem: filter_paired_n() first (the TMT ratio design
# leaves scattered missingness that would otherwise break the shared contrast matrix),
# technical covariates excluded, no voom weighting.
#
# Output: staging/freeze/proteomics/da/
#   human-precovid-sed-adu_{t07-adipose,t10-muscle}_prot-ph_da_dream-acute_v<x.y>.txt

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "09_build_da", "da_common.R"))

# ---- Upstream build, verbatim from .../generate_differential_analysis/generate_prot_ph_DA.R ----
#
# .generate_prot_ph_inputs = function(repo_local_dir,
#                                     model_type,
#                                     tissue,
#                                     parallel = F){
#   if(model_type == "training" & tissue == "adipose"){
#     message("Note: Adipose Prot-ph/pr has no training samples. No training analysis will be done.")
#   }else{
#     desired_ome = 'prot-ph'
#     ome_data = load_qc(selected_tissues = tissue,
#                        selected_omes = desired_ome,
#                        load_acute_only = FALSE)
#     ome_data = filter_paired_n(qc_data = ome_data,
#                                tissue = tissue,
#                                ome = desired_ome)
#
#     if (length(ome_data) > 0) {
#       metadata = ome_data[[tissue]][[desired_ome]][['sample_metadata']]
#       if (model_type == "acute") metadata = metadata %>% dplyr::filter(visitcode == 'ADU_BAS')  #filter just to the initial acute bout
#
#       rownames(metadata) = metadata$vialLabel
#       data_matrix = ome_data[[tissue]][[desired_ome]][['qc_norm']] %>%
#         dplyr::select(as.character(metadata$vialLabel))
#
#       process_metadata = process_covariates(meta = metadata,
#                                             selected_ome = desired_ome,
#                                             tissue_input = tissue,
#                                             include_technical = F)
#
#       .run_models(repo_local_dir = repo_local_dir,
#                   model_type = model_type,
#                   expression_object = data_matrix,
#                   process_metadata = process_metadata,
#                   tissue = tissue,
#                   ome = desired_ome,
#                   voom = FALSE,
#                   parallel = parallel)
#     }
#   }
# }
#
# Adapted below: load_qc() becomes load_qc_local() (this stage's step-08 *_QC.rda), and
# .run_models() becomes run_da_models() writing into staging/freeze. The training branch is
# dropped since only the acute analysis is built here. Modelling is unchanged.

# PARALLEL PATH: `parallel` reaches run_dream()'s MulticoreParam fork backend. See
# VARIANCEPARTITION_PARALLEL_CORES in config/pipeline.env for why it is fork and not SOCK.
generate_prot_ph_da = function(tissue, parallel = da_parallel_enabled()){
  desired_ome = "prot-ph"

  ome_data = load_qc_local(selected_tissues = tissue,
                           selected_omes = desired_ome,
                           load_acute_only = FALSE)
  if (length(ome_data) == 0)
    stop("no local *_QC object for ", tissue, "/", desired_ome, " — run step 08 first")

  ome_data = filter_paired_n(qc_data = ome_data, tissue = tissue, ome = desired_ome)

  metadata = ome_data[[tissue]][[desired_ome]][["sample_metadata"]]
  metadata = metadata %>% dplyr::filter(visitcode == "ADU_BAS")   # initial acute bout only
  rownames(metadata) = metadata$vialLabel

  data_matrix = ome_data[[tissue]][[desired_ome]][["qc_norm"]] %>%
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

for (tis in c("adipose", "muscle")) generate_prot_ph_da(tissue = tis)
