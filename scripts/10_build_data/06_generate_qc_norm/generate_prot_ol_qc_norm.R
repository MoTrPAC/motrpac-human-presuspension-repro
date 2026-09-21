#!/usr/bin/env Rscript
# Stage 1 QC-norm stem: prot-ol (Olink targeted proteomics)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_normalized_expression/
# generate_prot_ol_qc_norm.R.
#
# The in-house objects/functions it used from the installed packages — load_pheno,
# OUTLIERS, process_covariates, write_with_path_name — now come from the vendored
# lib/qc_helpers.R (which reads OME_TISSUE_CODE/OUTLIERS/COVARIATES_FILE from
# preflight and pheno from this stage's data/). The gated MotrpacBicQC::dl_read_gcp
# reads are unchanged. .annotate_olink hits biomaRt (Ensembl v‑online) for the
# feature-metadata output.
#
# Output: tab-delimited BIC-named files under
#   staging/freeze/proteomics/{metadata,qc-norm}/
# (files, not .rda — the *_QC objects are assembled from these freeze files via
# load_qc; that packaging step is separate.)

suppressWarnings(suppressMessages({ library(dplyr); library(tidyr); library(stringr) }))
# Root + shared internals via the PRECOVID_ROOT config env (config/pipeline.env).
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

generate_prot_ol_qc_norm = function(repo_local_dir){
  library(stringr)
  desired_ome = 'prot-ol'; tissue = 'blood'
  local_path = file.path(repo_local_dir)
  # Raw input paths resolved from the preflight quant-id catalog (not hard-coded).
  olink_data_path = quantid_path("prot-ol", tissue_code = "t02-plasma", data_details = "results")
  olink_data = MotrpacBicQC::dl_read_gcp(olink_data_path, tmpdir = .RAW_FILES, check_first = T)
  sample_metadata = MotrpacBicQC::dl_read_gcp(quantid_path("prot-ol", tissue_code = "t02-plasma", data_details = "metadata-samples"),
                                              tmpdir = .RAW_FILES, check_first = TRUE) %>%
    dplyr::rename(vialLabel = sample_id)
  protein_metadata = MotrpacBicQC::dl_read_gcp(quantid_path("prot-ol", tissue_code = "t02-plasma", data_details = "metadata-proteins"),
                                               tmpdir = .RAW_FILES, check_first = TRUE) %>%
    dplyr::rename(OlinkID = olink_id)

  metadata_path = paste0(local_path, "freeze/proteomics/metadata/")
  qc_norm_path = paste0(local_path, "freeze/proteomics/qc-norm/")
  dir.create(metadata_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_norm_path, recursive = TRUE, showWarnings = FALSE)

  #this filtering is borrowed from Jacob Barber's code
  message("Generating normalized matrixes for prot-ol")

  t_dat = olink_data %>% t() %>% as.data.frame()
  t_dat = t_dat %>% dplyr::mutate(vialLabel = rownames(t_dat)) %>% dplyr::select(vialLabel, everything()) #reorder vialLabel
  colnames(t_dat) = t_dat[1,] #set colnames then remove the first row
  t_dat = t_dat[-1,]

  t_dat = t_dat %>% dplyr::rename(vialLabel = olink_id) %>% dplyr::mutate_at(vars(starts_with('OID')), as.numeric)
  raw_olink = t_dat %>% tidyr::pivot_longer(cols=-vialLabel, names_to = 'OlinkID', values_to = 'NPX')
  raw_olink = raw_olink %>% dplyr::left_join(sample_metadata, by='vialLabel')
  raw_olink = raw_olink %>% dplyr::left_join(protein_metadata, by='OlinkID')
  raw_olink = raw_olink %>% dplyr::filter(missing_freq < 0.8)

  raw_olink_wide = raw_olink %>%
    dplyr::select("vialLabel","OlinkID","plate_id","NPX") %>%
    tidyr::pivot_wider(names_from = "OlinkID", values_from = "NPX") %>%
    na.omit()

  #here we remove outliers
  outliers_data = OUTLIERS$vialLabel
  pheno_data_parsed = load_pheno(load_acute_only = FALSE)$pheno_data

  raw_olink_wide_manifest <- dplyr::left_join(raw_olink_wide, pheno_data_parsed, by='vialLabel') %>%
    dplyr::select(pid, visitcode, plate_id, vialLabel, BMI, calculatedAge, Timepoint,study, sex_psca,starts_with('OID'))
  norm_adu <- raw_olink_wide_manifest %>%
    dplyr::filter(study=='01') %>%
    dplyr::filter(!vialLabel %in% as.character(outliers_data)) %>%
    dplyr::select(vialLabel, starts_with('OID'))

  new_norm_table <- norm_adu %>% dplyr::select(vialLabel, starts_with('OID')) %>% t()
  colnames(new_norm_table) <- new_norm_table[1,] #set viallabels as colnames
  new_norm_table <- new_norm_table[-1,] %>% as.data.frame() %>%
    dplyr::mutate(across(starts_with('1'), ~as.numeric(.)))

  sample_metadata_output = sample_metadata[sample_metadata$vialLabel %in% colnames(new_norm_table),]

  wider_metadata = merge(sample_metadata_output, pheno_data_parsed, by = 'vialLabel')
  # align metadata rows to the qc-norm matrix columns (samples) before batch
  # correction — merge() sorts by vialLabel, so without this removeBatchEffect
  # would apply each sample's covariates to the wrong column.
  wider_metadata = wider_metadata[match(colnames(new_norm_table), wider_metadata$vialLabel), ]
  process_metadata = process_covariates(meta = wider_metadata,
                                        selected_ome = desired_ome,
                                        tissue_input = tissue)
  meta = process_metadata$metadata
  technical_cov = paste(process_metadata[["technical_cov"]]$covariate, collapse = " + ")
  design_cov = paste(process_metadata[["design_cov"]], collapse = " + ")
  message(tissue, ";", desired_ome, ";technical: ", technical_cov, ";design: ", design_cov)

  new_norm_table = limma::removeBatchEffect(new_norm_table,
                                            covariates = model.matrix(as.formula(paste("~ ", technical_cov)), data = meta),
                                            design = model.matrix(as.formula(paste("~ ", design_cov)), data = meta))

  oids <- rownames(new_norm_table)
  new_norm_table <- as.data.frame(new_norm_table)
  new_norm_table$feature_id <- oids
  new_norm_table <- new_norm_table %>% dplyr::select(feature_id, everything()) #set feature_id as first col

  protein_metadata_output = .annotate_olink(protein_metadata) %>%
    dplyr::filter(feature_id %in% rownames(new_norm_table))

  write_with_path_name(sample_metadata_output, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'samples')
  write_with_path_name(new_norm_table, local_path = qc_norm_path, ome = desired_ome, tissue = tissue, data_category = 'qc-norm', data_details = 'log2')
  write_with_path_name(protein_metadata_output, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'features')
}

.annotate_olink = function(protein_metadata){
  #olink reshape; keep olink's gene target in "assay" as a check, but rename it generically
  olink = protein_metadata %>%
    dplyr::mutate(len_UP = stringr::str_count(uniprot_entry)) %>%
    tidyr::separate(uniprot_entry, into = c("uniprot", "redundant_ids"), sep = "_") %>%
    dplyr::mutate(gene_platform = stringr::str_remove(assay, pattern = "_.*")) %>%
    dplyr::select(feature_id = OlinkID, uniprot, gene_platform, redundant_ids) %>%
    dplyr::mutate(platform = "prot-ol")

  attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "uniprotswissprot")

  prot_lookup_df = ensembl_v105_getBM(attributes = attributes,
                                  filters = "uniprotswissprot",
                                  values = olink$uniprot) %>%
    dplyr::distinct() %>%
    dplyr::full_join(olink, by = c("uniprotswissprot" = "uniprot"))  %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(external_gene_name == "", NA, external_gene_name)) %>%
    dplyr::rename(ensembl_gene = ensembl_gene_id, entrez_gene = entrezgene_id, uniprot = uniprotswissprot) %>%
    dplyr::select(assay = platform, feature_id, entrez_gene, gene_symbol, ensembl_gene, uniprot) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::slice_min(entrez_gene, n = 1, with_ties = FALSE)

  #a few incorrect gene names, manually corrected upstream
  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "SHAN3"] <- "SHANK3"
  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "HECD4"] <- "HECTD4"
  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "WASH6"] <- "WASH6P"
  return(prot_lookup_df)
}

# Run: write the freeze files under the staging dir.
out_base <- .STAGING
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
generate_prot_ol_qc_norm(paste0(out_base, "/"))
message("prot-ol QC-norm written under ", file.path(out_base, "freeze", "proteomics"))
