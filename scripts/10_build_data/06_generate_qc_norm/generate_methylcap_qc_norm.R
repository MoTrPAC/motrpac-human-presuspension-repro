#!/usr/bin/env Rscript
# Stage 1 QC-norm stem: epigen-methylcap-seq (MethylCap-seq) — feature metadata only
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_normalized_expression/
# generate_methylcap_qc_norm.R.
#
# MethylCap-seq is NOT normalized in R — all preprocessing, QC and normalization
# run in an external, assay-specific pipeline produced by Yongchao Ge; that raw
# code is not included in this repo. This stem therefore writes no qc-norm matrix of its
# own: the vendored beta-value files ARE the qc-norm. What it does produce is the FEATURE
# METADATA (genomic-range -> gene annotation), consuming those matrices for their
# feature_ids (chr:start-end), and the SAMPLE METADATA built from the assay's
# qa-qc-metrics table.
#
# Those beta-value matrices are an external input, vendored under
# data-raw/sources/methylcap_qc_norm/ (one file per tissue) and read locally —
# see that folder's ../README.md for provenance. Their presence is enforced up
# front by scripts/10_build_data/check_required_inputs.sh, so this stem does no
# per-file existence checking of its own. In-house objects/functions
# (write_with_path_name, OME_TISSUE_CODE, the shared pre_cawg_get_peak_annotations_hs)
# come from the vendored lib/qc_helpers.R. Annotation uses biomaRt Ensembl v105 +
# ChIPseeker/txdbmaker (online).
#
# Output: BIC-named freeze files under staging/freeze/epigenomics/metadata/ —
#         metadata_features (version 1.4) and metadata_samples (default version, matching
#         how the ATAC stem versions its own sample metadata).
suppressWarnings(suppressMessages({ library(dplyr); library(tibble); library(data.table) }))
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

generate_methylcap_qc_norm = function(repo_local_dir){
  desired_ome = "epigen-methylcap-seq"
  metadata_path = file.path(repo_local_dir, "freeze", "epigenomics", "metadata")
  qc_norm_path = file.path(repo_local_dir, "freeze", "epigenomics", "qc-norm")
  dir.create(metadata_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_norm_path, recursive = TRUE, showWarnings = FALSE)

  # Assay-level sample metadata, exactly as the ATAC stem does it: one qa-qc-metrics table
  # covering every tissue, resolved from the preflight quant-id catalog and subset per tissue
  # below. Without this, methylcap was the only ome with no metadata_samples freeze, which
  # left step 08 with nothing to attach its beta-value matrices to.
  ome_meta = MotrpacBicQC::dl_read_gcp(quantid_path(desired_ome, data_details = "qa-qc-metrics"),
                                       tmpdir = .RAW_FILES, sep = ",", check_first = TRUE)
  ome_meta$vialLabel = as.character(ome_meta$vialLabel)
  pheno_data_parsed = load_pheno()$pheno_data

  # external artifact (Yongchao Ge; raw preprocessing not in this repo), vendored
  # locally and consumed only for feature_ids. Presence is guaranteed by the
  # Stage 1 required-inputs gate. See sources/methylcap_qc_norm/../README.md.
  src_dir = file.path(.QC_REPO, "scripts", "10_build_data", "06_generate_qc_norm", "sources", "methylcap_qc_norm")
  src_files = list.files(src_dir, pattern = "methylcap.*beta-values.*\\.txt$", full.names = TRUE)

  # tissue (blood/muscle/adipose) <-> tissue_code for methylcap, from the data object
  methyl_tissues = OME_TISSUE_CODE[OME_TISSUE_CODE$ome == desired_ome, c("tissue", "tissue_code")]

  for (file_path in src_files){
    matched = methyl_tissues$tissue_code[vapply(methyl_tissues$tissue_code,
                                                function(tc) grepl(tc, basename(file_path), fixed = TRUE),
                                                logical(1))]
    if (length(matched) != 1)
      stop("could not map ", basename(file_path), " to a single methylcap tissue_code")
    tissue = methyl_tissues$tissue[methyl_tissues$tissue_code == matched]
    message(paste("Generating feature metadata for", tissue))

    # There is no in-R normalization step here to write a matrix out of, so the vendored
    # beta-value file is copied into the freeze verbatim rather than round-tripped through
    # write_with_path_name(). Its filename already carries the standard BIC
    # <tissue_code>_<ome>_qc-norm token, so once copied it sits alongside every other ome's
    # qc-norm and the freeze tests, step 08 and everything downstream pick methylcap up
    # without needing a special case for it.
    file.copy(file_path, file.path(qc_norm_path, basename(file_path)), overwrite = TRUE)

    qc_norm = read.csv(file_path, sep = "\t", check.names = FALSE)
    # first column holds the feature_id (chr:start-end); the rest are vialLabels
    feature_metadata = data.frame(feature_id = qc_norm[[1]], stringsAsFactors = FALSE)
    methyl_annotation = .annotate_methyl_features(feature_metadata) %>%
      dplyr::arrange(feature_id)

    # restrict to the sed-adult cohort the same way ATAC does: pheno defines who is in
    # scope, and the qa-qc table is then subset to those vials.
    sample_cols = setdiff(colnames(qc_norm), colnames(qc_norm)[1])
    tissue_pheno = pheno_data_parsed[pheno_data_parsed$vialLabel %in% sample_cols, ]
    tissue_metadata = ome_meta[ome_meta$vialLabel %in% tissue_pheno$vialLabel, ]
    message(sprintf("  %s: %d beta-value sample(s), %d in pheno, %d with qa-qc metrics",
                    tissue, length(sample_cols), nrow(tissue_pheno), nrow(tissue_metadata)))

    write_with_path_name(
      actual_data_object = tissue_metadata,
      local_path = metadata_path,
      tissue = tissue,
      ome = desired_ome,
      data_category = "metadata",
      data_details = "samples"
    )

    write_with_path_name(
      actual_data_object = methyl_annotation,
      local_path = metadata_path,
      tissue = tissue,
      ome = desired_ome,
      data_category = "metadata",
      data_details = "features",
      return_name_only = FALSE
    )
  }
}

.annotate_methyl_features = function(feature_metadata){
  feature_input = feature_metadata %>%
    dplyr::mutate(chrom = gsub(":.*", "", feature_id),
                  start = as.numeric(gsub(".*:|-.*", "", feature_id)),
                  end = as.numeric(gsub(".*-", "", feature_id))) %>%
    dplyr::mutate(chrom = gsub("chrM", "chrMT", chrom)) %>%
    data.table::as.data.table()

  methyl_peakdf = pre_cawg_get_peak_annotations_hs(feature_input)
  attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name")
  methyl_lookup_df = ensembl_v105_getBM(attributes = attributes,
                                    filters = "ensembl_gene_id",
                                    values = methyl_peakdf$ensembl_gene) %>%
    dplyr::full_join(methyl_peakdf, c("ensembl_gene_id" = "ensembl_gene")) %>%
    dplyr::mutate(gene_symbol = dplyr::na_if(external_gene_name, "")) %>%
    dplyr::select(feature_id,
                  entrez_gene = entrezgene_id,
                  gene_symbol,
                  ensembl_gene = ensembl_gene_id,
                  custom_annotation,
                  relationship_to_gene) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::slice_min(entrez_gene, n = 1, with_ties = FALSE) %>%
    dplyr::mutate(entrez_gene = as.character(entrez_gene),
                  assay = "epigen-methylcap-seq") %>%
    dplyr::relocate(assay)

  return(methyl_lookup_df)
}

out_base <- .STAGING
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
generate_methylcap_qc_norm(out_base)
message("MethylCap-seq feature metadata written under ", file.path(out_base, "freeze", "epigenomics", "metadata"))
