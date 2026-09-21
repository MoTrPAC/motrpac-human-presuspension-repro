#!/usr/bin/env Rscript
# Stage 1 QC-norm stem: transcript-rna-seq (RNA-seq)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_normalized_expression/
# generate_transcriptomics_qc_norm.R.
#
# In-house objects/functions (load_pheno, OUTLIERS, process_covariates,
# write_with_path_name) come from the vendored lib/qc_helpers.R. Raw qa-qc-metrics
# and rsem count paths are resolved from the preflight quant-id catalog via
# quantid_path() (the upstream hard-coded human-precovid/results paths are dead).
# Gated dl_read_gcp + edgeR/limma + biomaRt (Ensembl v105, online, feature
# metadata) are unchanged. Adds the metadata<->matrix match() alignment before
# batch correction (upstream fed a vialLabel-sorted meta against count-file-order
# columns — the same misalignment fixed in prot-ol).
#
# Output: BIC-named freeze files under staging/freeze/transcriptomics/.
suppressWarnings(suppressMessages({ library(dplyr); library(tibble); library(stringr) }))
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

generate_transcriptomics_qc_norm = function(repo_local_dir){
  require("biomaRt")
  desired_ome = 'transcript-rna-seq'
  tissue_types = c('muscle', 'blood', 'adipose')
  local_path = repo_local_dir

  merged_metadata = MotrpacBicQC::dl_read_gcp(quantid_path(desired_ome, data_details = "qa-qc-metrics"),
                                              tmpdir = .RAW_FILES, sep = ",", check_first = TRUE) %>%
    dplyr::select(-"BID") %>%
    dplyr::select(-"PID")

  merged_metadata$Batch <- sub(".*?(\\d{1,2})$", "\\1", merged_metadata$Lib_batch_ID) #choose last 2 characters
  pheno_data_parsed = load_pheno(load_acute_only = FALSE)$pheno_data

  metadata_path = paste0(local_path, "freeze/transcriptomics/metadata/")
  qc_norm_path = paste0(local_path, "freeze/transcriptomics/qc-norm/")
  dir.create(metadata_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_norm_path, recursive = TRUE, showWarnings = FALSE)

  for (tissue in tissue_types){
    message(paste("Generating transcriptomics normalized matrixes for", tissue))
    tc = .match_ome_tissue_code(desired_ome, tissue)
    file_load = quantid_path(desired_ome, tissue_code = tc, data_details = "rsem-genes-count")
    raw_counts_input = MotrpacBicQC::dl_read_gcp(file_load, sep = '\t', tmpdir = .RAW_FILES) %>%
      tibble::column_to_rownames("gene_id")

    outliers_vialLabels = OUTLIERS$vialLabel
    desired_samples = setdiff(colnames(raw_counts_input), outliers_vialLabels) #counts minus outliers
    tissue_pheno = pheno_data_parsed[pheno_data_parsed$vialLabel %in% desired_samples, ]
    tissue_metadata = merged_metadata[merged_metadata$vialLabel %in% tissue_pheno$vialLabel, ]
    raw_counts_input = raw_counts_input[, colnames(raw_counts_input) %in% tissue_pheno$vialLabel] #remove HA, peds

    meta = merge(tissue_pheno, tissue_metadata, by = "vialLabel")
    # align counts columns and metadata rows before batch correction: keep only
    # metadata-covered columns, then order meta to the count-matrix columns.
    raw_counts_input = raw_counts_input[, colnames(raw_counts_input) %in% meta$vialLabel]
    meta = meta[match(colnames(raw_counts_input), meta$vialLabel), ]
    process_metadata = process_covariates(meta = meta,
                                          selected_ome = desired_ome,
                                          tissue_input = tissue)
    meta = process_metadata$metadata

    #-----filter lowly expressed genes and transform into logcpm
    raw_dge = edgeR::DGEList(counts = raw_counts_input)
    keep = rowSums(edgeR::cpm(raw_dge) > 0.5) >= round(length(raw_counts_input)*0.1)
    filt_dge = raw_dge[keep, , keep.lib.sizes=FALSE]
    dge = edgeR::calcNormFactors(filt_dge, method="TMM")
    norm_counts = edgeR::cpm(dge,log=TRUE)

    technical_cov = paste(process_metadata[["technical_cov"]]$covariate, collapse = " + ")
    design_cov = paste(process_metadata[["design_cov"]], collapse = " + ")
    message(tissue," technical: ", technical_cov, " design: ", design_cov)
    #---batch correction :: ONLY for visualization/clustering/etc. Use raw counts for DA
    batch_corrected = limma::removeBatchEffect(norm_counts,
                                               covariates = stats::model.matrix(stats::as.formula(paste("~ ", technical_cov)), data = meta),
                                               design = stats::model.matrix(stats::as.formula(paste("~ ", design_cov)), data = meta))
    batch_corrected = as.data.frame(batch_corrected)
    batch_corrected$feature_id = rownames(batch_corrected)
    batch_corrected = batch_corrected %>% dplyr::select(feature_id, everything())

    only_features = batch_corrected %>% dplyr::select(feature_id)
    rna_gene_anno = .annotate_rna_features(only_features)

    write_with_path_name(tissue_metadata, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'samples')
    write_with_path_name(batch_corrected, local_path = qc_norm_path, ome = desired_ome, tissue = tissue, data_category = 'qc-norm', data_details = 'log-cpm')
    #version 1.4 adds the gene symbols directly to the metadata-features.
    write_with_path_name(rna_gene_anno, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'features')
  }
}

# Requires dbplyr 2.3.4 (BioMart lazy-table collect issue). Ensembl v105 / GENCODE 39.
.annotate_rna_features = function(only_features){
  rna_features <- only_features %>%
    dplyr::mutate(ensembl_gene_id = stringr::str_remove(feature_id, "\\..*")) %>%
    dplyr::distinct()
  attributes <- c("ensembl_gene_id", "entrezgene_id", "external_gene_name")
  rna_lookup_df <- ensembl_v105_getBM(attributes = attributes,
                                  filters = "ensembl_gene_id",
                                  values = rna_features$ensembl_gene_id) %>%
    dplyr::full_join(rna_features, by = "ensembl_gene_id", relationship = "many-to-many") %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(external_gene_name == "", NA, external_gene_name)) %>%
    dplyr::select(feature_id,
                  entrez_gene = entrezgene_id,
                  gene_symbol,
                  ensembl_gene = ensembl_gene_id) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::slice_min(entrez_gene, n = 1, with_ties = FALSE) %>%
    dplyr::mutate(entrez_gene = as.character(entrez_gene),
                  assay = "transcript-rna-seq") %>%
    dplyr::relocate(assay) %>%
    dplyr::distinct()
  return(rna_lookup_df)
}

out_base <- .STAGING
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
generate_transcriptomics_qc_norm(paste0(out_base, "/"))
message("transcriptomics QC-norm written under ", file.path(out_base, "freeze", "transcriptomics"))
