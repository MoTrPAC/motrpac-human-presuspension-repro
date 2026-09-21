#!/usr/bin/env Rscript
# Step 09 DA stem: metabolomics (all targeted + untargeted platforms, all three tissues).
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_differential_analysis/
# generate_metabolomics_DA.R.
#
# One DA table per tissue x PLATFORM — the platforms are pooled into a single
# {TISSUE}_METAB_DA object only later, at the step-10 assembly. Technical covariates are
# included here (unlike the proteomics stems) and no voom weighting is applied.
#
# load_qc_local() is called once PER PLATFORM, not once per tissue with every platform
# selected. That is deliberate: remove_redundant_metab (on by default) runs
# .eliminate_redundant_metab(), which drops duplicated metabolites by lowest CV ACROSS
# whatever is loaded together. Loading all platforms in one call would dedupe between
# platforms and silently change the feature set, so the per-platform loop matches upstream.
#
# Output: staging/freeze/metabolomics-{targeted,untargeted}/da/
#   human-precovid-sed-adu_<tissue_code>_<platform>_da_dream-acute_v<x.y>.txt

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "09_build_da", "da_common.R"))

# ---- Upstream build, verbatim from .../generate_differential_analysis/generate_metabolomics_DA.R ----
#
# .generate_metabolomics_inputs = function(repo_local_dir,
#                                          model_type,
#                                          tissue,
#                                          assay,
#                                          parallel = F){
#   message(assay); message(tissue)
#   desired_ome = assay
#   ome_data = load_qc(selected_tissue = tissue,
#                      selected_ome = desired_ome,
#                      load_acute_only = FALSE,
#                      remove_unnamed_metab = TRUE)
#   if (length(ome_data) > 0 && nrow(ome_data[[tissue]][[desired_ome]][["qc_norm"]] > 0)) {
#     metadata = ome_data[[tissue]][[desired_ome]][['sample_metadata']]
#     if (model_type == "acute") metadata = metadata %>% dplyr::filter(visitcode == 'ADU_BAS')  #filter just to the initial acute bout
#     rownames(metadata) = metadata$vialLabel
#     data_matrix = ome_data[[tissue]][[desired_ome]][['qc_norm']] %>%
#       dplyr::select(as.character(metadata$vialLabel))
#
#     process_metadata = process_covariates(meta = metadata,
#                                           selected_ome = desired_ome,
#                                           tissue_input = tissue,
#                                           include_technical = TRUE)
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
# Adapted below: load_qc() becomes load_qc_local() (this stage's step-08 *_QC.rda) with the
# plural argument names spelled out — upstream relies on R partial matching of
# selected_tissue/selected_ome — and .run_models() becomes run_da_models() writing into
# staging/freeze. The upstream emptiness guard reads `nrow(x > 0)`, which parenthesises the
# comparison inside nrow() and happens to work; it is written explicitly here. Modelling is
# unchanged.
#
# Available metab platforms for a tissue, derived from the local *_QC objects using the same
# object-name -> (tissue, ome) mapping as load_qc_local(). Clinical chemistry has its own
# upstream per-feature DA generator, so it is excluded here.
.available_metab_omes = function(tissue){
  f = list.files(file.path(ROOT, "scripts", "10_build_data", "data"), pattern = "_QC\\.rda$")
  obj = sub("\\.rda$", "", f)
  tis = tolower(sub("_.*", "", obj))
  omes = gsub("_", "-", tolower(sub("^[^_]+_(.*)_QC$", "\\1", obj)))
  out = omes[tis == tissue & grepl("^metab-", omes)]
  return(sort(setdiff(out, "metab-t-clinical")))
}

# PARALLEL PATH: `parallel` reaches run_dream()'s MulticoreParam fork backend. See
# VARIANCEPARTITION_PARALLEL_CORES in config/pipeline.env for why it is fork and not SOCK.
generate_metab_da = function(tissue, desired_ome, parallel = da_parallel_enabled()){
  message(desired_ome); message(tissue)

  # remove_redundant_metab = FALSE, overriding load_qc_local()'s default. See the adaptation
  # note above: the CV-based dedup belongs downstream of the DA, not in front of it.
  ome_data = load_qc_local(selected_tissues = tissue,
                           selected_omes = desired_ome,
                           load_acute_only = FALSE,
                           remove_unnamed_metab = TRUE,
                           remove_redundant_metab = FALSE)
  if (length(ome_data) == 0) {
    warning("no local *_QC object for ", tissue, "/", desired_ome, " — skipping")
    return(invisible(NULL))
  }
  leaf = ome_data[[tissue]][[desired_ome]]
  if (is.null(leaf) || nrow(leaf[["qc_norm"]]) == 0) {
    warning("empty qc_norm for ", tissue, "/", desired_ome, " — skipping")
    return(invisible(NULL))
  }

  metadata = leaf[["sample_metadata"]]
  metadata = metadata %>% dplyr::filter(visitcode == "ADU_BAS")   # initial acute bout only
  rownames(metadata) = metadata$vialLabel

  data_matrix = leaf[["qc_norm"]] %>% dplyr::select(as.character(metadata$vialLabel))

  process_metadata = process_covariates(meta = metadata,
                                        selected_ome = desired_ome,
                                        tissue_input = tissue,
                                        include_technical = TRUE)

  # run_dream() rather than run_da_models(): the latter writes immediately, and the
  # platform/assay columns below have to be added before the file is written.
  # Sample alignment is explicit inside run_dream(), which reorders the covariate rows with
  # match(colnames(expression_object), rownames(meta_matrix)) before fitting, so the design
  # always follows the expression columns rather than their incoming order.
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

  # .convert_dream_output() sets assay = ome for every stem. For metab the ome IS the
  # platform, so carry both: assay = "metab" (the value HUMAN_FEATURE_TO_GENE keys on) and
  # platform = the ome. Adapted from differential_analysis_results.R:73-80, minus the
  # pooling — the freeze keeps one file per tissue x platform, as the release does.
  out = out %>%
    dplyr::mutate(platform = assay, assay = "metab") %>%
    dplyr::relocate(platform, .after = assay)

  da_path = .da_freeze_dir(desired_ome)
  write_with_path_name(out,
                       local_path = da_path,
                       ome = desired_ome,
                       tissue = tissue,
                       data_category = "da",
                       data_details = paste0("dream-", DA_MODEL_TYPE))
  message(sprintf("%s / %s: %d rows across %d contrast(s) -> %s",
                  tissue, desired_ome, nrow(out), dplyr::n_distinct(out$contrast), da_path))
  return(invisible(out))
}

for (tis in c("adipose", "blood", "muscle")) {
  platforms = .available_metab_omes(tis)
  message(sprintf("== %s: %d metab platform(s): %s ==", tis, length(platforms),
                  paste(platforms, collapse = ", ")))
  for (plat in platforms) generate_metab_da(tissue = tis, desired_ome = plat)
}
