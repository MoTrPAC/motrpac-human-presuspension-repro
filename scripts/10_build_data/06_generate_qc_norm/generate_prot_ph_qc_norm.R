#!/usr/bin/env Rscript
# Stage 1 QC-norm stem: prot-ph (phosphoproteomics)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_normalized_expression/
# generate_prot_ph_qc_norm.R.
#
# In-house objects/functions (load_pheno, OUTLIERS, process_covariates,
# write_with_path_name, .median_mad_norm, .remove_na, report_replicate_correlations)
# come from the vendored lib/qc_helpers.R. Raw ratio-results + vial-metadata paths
# are resolved from the preflight quant-id catalog via quantid_path() (the upstream hard-coded
# human-precovid/results paths are dead). Gated dl_read_gcp + cmapR (GCT) + biomaRt
# (online, feature metadata) are unchanged. This stem already aligns metadata to
# the matrix columns via match()/GCT before batch correction.
#
# Output: BIC-named freeze files under staging/freeze/proteomics/.
suppressWarnings(suppressMessages({ library(dplyr); library(stringr) }))
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

generate_prot_ph_qc_norm = function(repo_local_dir){
  local_path = repo_local_dir
  desired_ome = 'prot-ph'; tissue_types = c('muscle', 'adipose')

  pheno_data_parsed = load_pheno(load_acute_only = FALSE)$pheno_data
  metadata_path = paste0(local_path, "freeze/proteomics/metadata/")
  qc_norm_path = paste0(local_path, "freeze/proteomics/qc-norm/")
  dir.create(metadata_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_norm_path, recursive = TRUE, showWarnings = FALSE)

  for (tissue in tissue_types){
    message(paste("Generating normalized matrixes for phospho-ph", tissue))
    tc = .match_ome_tissue_code(desired_ome, tissue)
    file_load = quantid_path(desired_ome, tissue_code = tc, data_details = "ratio-results")
    phospho = MotrpacBicQC::dl_read_gcp(file_load, tmpdir = .RAW_FILES)
    phospho = phospho[!duplicated(phospho$ptm_id),] #remove duplicates
    phospho = phospho %>% dplyr::filter(phospho$is_contaminant == "FALSE")

    mat = phospho %>% dplyr::select(matches("^[0-9]")) %>% as.data.frame()
    rdesc <- phospho %>% dplyr::select(!matches("^[0-9]"))
    rdesc <- rdesc %>% dplyr::select(!starts_with("pool"))

    rdesc$id <- rdesc$ptm_id
    rownames(mat) = rdesc$ptm_id

    tmt_metadata = MotrpacBicQC::dl_read_gcp(quantid_path(desired_ome, tissue_code = tc, data_details = "vial-metadata"), tmpdir = .RAW_FILES) %>%
      dplyr::rename(vialLabel= vial_label) %>%
      dplyr::mutate(id = vialLabel) %>%
      dplyr::mutate(vialLabel = gsub("\\.1", "", vialLabel))
    meta_merge = dplyr::left_join(tmt_metadata, pheno_data_parsed, by = 'vialLabel')
    meta_merge <- meta_merge[match(colnames(mat),meta_merge$id),]

    if(ncol(mat)!=nrow(meta_merge) | nrow(mat)!=nrow(rdesc)){
      stop("Column or row annotations do not match
		     the matrix size.")}

    #this is non-batch effect corrected, non-normalized dataset
    phospho_nonnorm <- cmapR::GCT(mat=as.matrix(mat), rdesc=rdesc, cdesc=meta_merge,
                                  rid =rownames(mat), cid = colnames(mat))
    phospho_nonnorm <- cmapR::subset_gct(phospho_nonnorm, cid = which(phospho_nonnorm@cdesc$study == "01")) #remove HA participants
    #This is non-batch effect corrected, median-normalized dataset
    phospho_mednorm <- .median_mad_norm(phospho_nonnorm, mad = F)

    remove = OUTLIERS$vialLabel
    phospho_mednorm_no_outliers <- cmapR::subset_gct(phospho_mednorm, cid=which(!(phospho_mednorm@cdesc$vialLabel %in% remove)))

    #batch correction
    phospho_sample_meta = phospho_mednorm_no_outliers@cdesc
    process_metadata = process_covariates(meta = phospho_sample_meta,
                                          selected_ome = desired_ome,
                                          tissue_input = tissue)
    meta = process_metadata$metadata
    technical_cov = paste(process_metadata[["technical_cov"]]$covariate, collapse = " + ")
    design_cov = paste(process_metadata[["design_cov"]], collapse = " + ")
    message(tissue, ";", desired_ome, ";technical: ", technical_cov, ";design: ", design_cov)

    phospho_mednorm_no_outliers@mat = limma::removeBatchEffect(phospho_mednorm_no_outliers@mat,
                                                               covariates = model.matrix(as.formula(paste("~ ", technical_cov)), data = meta),
                                                               design = model.matrix(as.formula(paste("~ ", design_cov)), data = meta))
    if(tissue == 'muscle'){
      matr <- phospho_mednorm_no_outliers@mat %>% as.data.frame(check.names = FALSE)
      # NOT IMPLEMENTED CORRECTLY BEFORE v2.0. Releases up to and including v1.3 wrote the
      # intra-site average into the ".1" column and then deleted that column, so the surviving
      # column kept its first measurement alone: the second measurement was dropped rather
      # than merged. Inter-site pairs were always averaged correctly, which is what made it
      # easy to miss. v2.0 is the first release in which intra-site replicates are merged.
      # 16 sample pairs are affected in muscle prot-ph and 16 in muscle prot-pr.
      #
      # Muscle carries two kinds of repeat measurement of the same sample: INTRA-site, where
      # the repeat's id is the original plus a ".1" suffix, and INTER-site, where the repeat's
      # id ends in 7 against a partner ending in 2. Each pair is merged by averaging, then the
      # redundant column is dropped.
      #
      # The mean is written into the column that SURVIVES the drop, which for an intra-site
      # pair is the un-suffixed id. Writing it into the ".1" column instead would leave the
      # surviving column holding the first measurement alone and discard the second — the
      # average would be computed and then deleted with the column it was written to.
      #
      # `\\.1$` is escaped: an unescaped ".1$" also matches any id ending in a digit then 1,
      # pulling non-replicate samples into the list.
      #
      # Each pair is derived once, into `kept` (the surviving column, which receives the mean)
      # and `dropped`, so the correlation report and the averaging can never disagree about
      # which two columns form a pair.
      intrasite_pairs = data.frame(dropped = phospho_mednorm_no_outliers@cdesc$id[grepl("\\.1$", phospho_mednorm_no_outliers@cdesc$id)])
      intrasite_pairs$kept = gsub("\\.1$", "", intrasite_pairs$dropped)
      #replicates across chemical sites
      intersite_pairs = data.frame(dropped = phospho_mednorm_no_outliers@cdesc$id[grepl("7$", phospho_mednorm_no_outliers@cdesc$id)])
      intersite_pairs$kept = paste0(stringr::str_sub(intersite_pairs$dropped, end = -2), "2")

      # Reported before either loop runs, so every correlation is between two raw measurements.
      report_replicate_correlations(matr, intrasite_pairs, pair_type = "intra-site",
                                    ome = desired_ome, tissue = tissue)
      report_replicate_correlations(matr, intersite_pairs, pair_type = "inter-site",
                                    ome = desired_ome, tissue = tissue)

      for(k in seq_len(nrow(intrasite_pairs))){
        matr[intrasite_pairs$kept[k]] <- rowMeans(matr[, c(intrasite_pairs$dropped[k], intrasite_pairs$kept[k])], na.rm = TRUE)
      }
      for(k in seq_len(nrow(intersite_pairs))){
        matr[intersite_pairs$kept[k]] <- rowMeans(matr[, c(intersite_pairs$kept[k], intersite_pairs$dropped[k])], na.rm = TRUE)
      }
      #removing columns which have been averaged
      matr <- dplyr::select(matr, -all_of(intrasite_pairs$dropped))
      matr <- dplyr::select(matr, -all_of(intersite_pairs$dropped))
      matr <- as.matrix(matr)

      cdescc = phospho_sample_meta %>%
        dplyr::filter(!grepl("\\.1$", id)) %>%
        dplyr::filter(id %in% colnames(matr))
      cdescc <- cdescc[match(colnames(matr),cdescc$id),]

      if(ncol(matr)!=nrow(cdescc)|nrow(mat)!=nrow(rdesc)){
        stop("Column or row annotations do not match
		     the matrix size. (after avg)")}

      phospho_ave <- cmapR::GCT(mat=matr, rdesc= data.frame(rdesc), cdesc=data.frame(cdescc), rid =rownames(matr), cid = colnames(matr))
      PH <- phospho_ave %>% .remove_na(0.7)
    }else{
      PH = phospho_mednorm_no_outliers %>% .remove_na(0.7)
    }
    phospho_output <- as.data.frame(PH@mat)
    phospho_output$feature_id = rownames(phospho_output)
    phospho_output = phospho_output %>% dplyr::select(feature_id, everything())

    sample_metadata_output = PH@cdesc %>% dplyr::select(c("vialLabel", "tmt_plex", "tmt16_channel"))

    feature_metadata_output = PH@rdesc %>%
      dplyr::rename(feature_id = ptm_id) %>%
      dplyr::select(feature_id, everything()) %>%
      dplyr::filter(feature_id %in% rownames(phospho_output)) %>%
      .annotate_prot_ph()

    write_with_path_name(sample_metadata_output, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'samples')
    write_with_path_name(phospho_output, local_path = qc_norm_path, ome = desired_ome, tissue = tissue, data_category = 'qc-norm', data_details = 'log2-mn')
    write_with_path_name(feature_metadata_output, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'features')
  }
}

.annotate_prot_ph = function(feature_metadata_output){
  #feature_id is the ptm_id; protein_id is the UniProt accession. Strip isoform suffixes for BioMart.
  prot_ph = feature_metadata_output %>%
    dplyr::mutate(platform = "prot-ph", uniprot = protein_id,
                  uniprot_lookup = stringr::str_remove(uniprot, "-.*"))

  attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "uniprotswissprot")

  prot_lookup_df = ensembl_v105_getBM(attributes = attributes, filters = "uniprotswissprot",
                                  values = prot_ph$uniprot_lookup) %>%
    dplyr::distinct() %>%
    dplyr::full_join(prot_ph, by = c("uniprotswissprot" = "uniprot_lookup")) %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(external_gene_name == "", NA, external_gene_name)) %>%
    dplyr::rename(ensembl_gene = ensembl_gene_id, entrez_gene = entrezgene_id) %>%
    # confident_site: the PTM site-localization flag, carried straight through from the raw
    # ratio-results file (it is already on `rdesc`, so this select() is the only thing that
    # ever removed it). Kept because it is the filter PTM-SEA input is built on — step 15
    # emits only fully localized sites, and an ambiguous site scored as if it were localized
    # attributes a phosphorylation to the wrong residue. It belongs beside flanking_sequence
    # for the same reason that one is kept: prot-ph-only feature annotation that a downstream
    # step cannot reconstruct from the matrix.
    dplyr::select(assay = platform, feature_id, entrez_gene, gene_symbol, ensembl_gene, uniprot, flanking_sequence, confident_site) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::slice_min(entrez_gene, n = 1, with_ties = FALSE)

  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "SHAN3"] <- "SHANK3"
  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "HECD4"] <- "HECTD4"
  prot_lookup_df$gene_symbol[prot_lookup_df$gene_symbol == "WASH6"] <- "WASH6P"
  return(prot_lookup_df)
}

out_base <- .STAGING
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
generate_prot_ph_qc_norm(paste0(out_base, "/"))
message("prot-ph QC-norm written under ", file.path(out_base, "freeze", "proteomics"))
