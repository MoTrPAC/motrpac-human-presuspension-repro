#!/usr/bin/env Rscript
# Step 09 DA stem: transcript-rna-seq (blood, muscle, adipose).
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_differential_analysis/
# generate_transcriptomics_DA.R.
#
# The only stem in this scope that does NOT model the qc-norm matrix. RNA-seq is fit on RAW
# COUNTS: the qc-norm object supplies the gene and sample selection, then an edgeR DGEList is
# built from the raw rsem-genes-count file and passed to dream with voom weighting. Technical
# covariates ARE included here (unlike the proteomics stems), matching upstream.
#
# Raw counts are resolved from the preflight quant-id catalog via quantid_path() rather than
# the upstream live `gsutil ls` (.find_path_name was dropped from qc_helpers for that reason),
# and cached under staging/raw-files by dl_read_gcp(check_first = TRUE).
#
# Output: staging/freeze/transcriptomics/da/
#   human-precovid-sed-adu_{t04-blood-rna,t06-muscle,t11-adipose}_transcript-rna-seq_da_dream-acute_v<x.y>.txt

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "09_build_da", "da_common.R"))
suppressWarnings(suppressMessages({library(tibble)}))

# ---- Upstream build, verbatim from .../generate_differential_analysis/generate_transcriptomics_DA.R ----
#
# .generate_transcriptomics_inputs = function(repo_local_dir,
#                                             model_type,
#                                             tissue,
#                                             parallel = F){
#   # check_package_installation("edgeR")
#   desired_ome = 'transcript-rna-seq'
#   local_path = paste0(repo_local_dir, "data/tmp/")
#   counts_data_path = .find_path_name(desired_ome = desired_ome, tissue = tissue, data_type = "rsem-genes-count")
#   if(length(counts_data_path) == 0) stop(paste("The desired", tissue, "is not available for transcriptomics"))
#
#   raw_counts_input = MotrpacBicQC::dl_read_gcp(counts_data_path, sep = '\t', tmpdir = local_path)
#   parsed_qc_norm = load_qc(selected_omes = desired_ome,
#                            selected_tissues = tissue,
#                            load_acute_only = FALSE)
#   #-> so here we want to make sure the genes we chose and the participants we chose are the same as the ones that pass through expression logcpm cutoffs
#   metadata = parsed_qc_norm[[tissue]][[desired_ome]][['sample_metadata']]
#   if (model_type == "acute") metadata = metadata %>% dplyr::filter(visitcode == 'ADU_BAS')  #filter just to the initial acute bout
#   rownames(metadata) = metadata$vialLabel
#
#   raw_counts_input = raw_counts_input %>%
#     dplyr::filter(gene_id %in% rownames(parsed_qc_norm[[tissue]][[desired_ome]][['qc_norm']]))  %>% #filter to select genes
#     tibble::column_to_rownames("gene_id") %>% #set to rownames
#     dplyr::select(as.character(metadata$vialLabel)) #reorganize and filter to only desired participants
#
#   pre_dge <- edgeR::DGEList(counts = raw_counts_input) #standard dge processing
#   pre_dge <- edgeR::calcNormFactors(pre_dge) #standard dge processing
#   process_metadata = process_covariates(meta = metadata,
#                                         selected_ome = desired_ome,
#                                         tissue_input = tissue)
#
#   .run_models(repo_local_dir = repo_local_dir,
#               model_type = model_type,
#               expression_object = pre_dge,
#               process_metadata = process_metadata,
#               tissue = tissue,
#               ome = desired_ome,
#               voom = TRUE,
#               parallel = parallel)
#
# }
#
# Adapted below: load_qc() becomes load_qc_local(), .find_path_name()'s live gsutil listing
# becomes quantid_path() over the preflight catalog, and .run_models() becomes
# run_da_models() writing into staging/freeze. Modelling is unchanged.

# PARALLEL PATH: `parallel` reaches run_dream()'s MulticoreParam fork backend. See
# VARIANCEPARTITION_PARALLEL_CORES in config/pipeline.env for why it is fork and not SOCK.
generate_transcriptomics_da = function(tissue, parallel = da_parallel_enabled()){
  desired_ome = "transcript-rna-seq"
  tissue_code = .match_ome_tissue_code(desired_ome = desired_ome, input_tissue = tissue)

  counts_data_path = quantid_path(desired_ome, tissue_code = tissue_code,
                                  data_details = "rsem-genes-count")
  if (length(counts_data_path) == 0)
    stop("no rsem-genes-count path in the quant-id catalog for ", tissue, "/", desired_ome)

  raw_counts_input = MotrpacBicQC::dl_read_gcp(counts_data_path, sep = "\t",
                                               tmpdir = .RAW_FILES, check_first = TRUE)

  parsed_qc_norm = load_qc_local(selected_omes = desired_ome,
                                 selected_tissues = tissue,
                                 load_acute_only = FALSE)
  if (length(parsed_qc_norm) == 0)
    stop("no local *_QC object for ", tissue, "/", desired_ome, " — run step 08 first")

  metadata = parsed_qc_norm[[tissue]][[desired_ome]][["sample_metadata"]]
  metadata = metadata %>% dplyr::filter(visitcode == "ADU_BAS")   # initial acute bout only
  rownames(metadata) = metadata$vialLabel

  # the qc-norm object defines which genes and samples survived the log-CPM cutoffs;
  # the model is then fit on the raw counts for exactly that selection.
  raw_counts_input = raw_counts_input %>%
    dplyr::filter(gene_id %in% rownames(parsed_qc_norm[[tissue]][[desired_ome]][["qc_norm"]])) %>%
    tibble::column_to_rownames("gene_id") %>%
    dplyr::select(as.character(metadata$vialLabel))

  pre_dge = edgeR::DGEList(counts = raw_counts_input)
  pre_dge = edgeR::calcNormFactors(pre_dge)

  process_metadata = process_covariates(meta = metadata,
                                        selected_ome = desired_ome,
                                        tissue_input = tissue)

  # Sample alignment is explicit inside run_da_models(): run_dream() reorders the covariate
  # rows with match(colnames(expression_object), rownames(meta_matrix)) before fitting, so
  # the design always follows the expression columns rather than their incoming order.
  out = run_da_models(expression_object = pre_dge,
                      process_metadata = process_metadata,
                      tissue = tissue,
                      ome = desired_ome,
                      voom = TRUE,
                      parallel = parallel)
  return(invisible(out))
}

for (tis in c("blood", "muscle", "adipose")) generate_transcriptomics_da(tissue = tis)
