#!/usr/bin/env Rscript
# Stage 1 QC-norm stem: clinical chemistry (prot-clinical + metab-t-clinical)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_normalized_expression/
# generate_clinical_qc_norm.R (author: Christopher Jin).
#
# Clinical chemistry is NOT normalized here — the analyte values are curated and
# normalized upstream by the MoTrPAC clinical group and ingested as-is (absolute
# terms, no log2). This stem's real work is the FEATURE METADATA: the 9 analytes
# split across two omes by molecule class (prot-clinical -> proteomics/,
# metab-t-clinical -> metabolomics-targeted/), annotated and written into the
# freeze layout alongside the qc-norm (absolute) and sample-metadata files.
#
# In-house objects/functions come from local .rda + the vendored lib/qc_helpers.R:
#   - clinical chemistry tables: the local cln_chemistry_* data objects in this
#     stage's data/ (were MotrpacHumanPreSuspensionData::load_clinical_data()[["chemistry"]])
#   - sample metadata: local pheno$data (was MotrpacHumanPreSuspensionData::pheno$data)
#   - write_with_path_name / .clinical_ome_for_feature / .freeze_subdir_for_ome
# Annotation is online: biomaRt (prot analytes -> gene) + the RefMet API (metab
# analytes -> refmet/kegg).
#
# Output: BIC-named freeze files under staging/freeze/{proteomics,metabolomics-targeted}/
#         (qc-norm + feature/sample metadata, version 1.4).
suppressWarnings(suppressMessages({ library(dplyr); library(tibble); library(stringr) }))
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

generate_clinical_qc_norm = function(repo_local_dir){
  freeze_output = file.path(repo_local_dir, "freeze")

  # clinical chemistry tables from the local cln_chemistry_* data objects
  # (were load_clinical_data()[["chemistry"]]). Each is a data.frame keyed by
  # analyte_name with per-sample (vialLabel) columns.
  chem_files = list.files(.STAGE1_DATA, pattern = "^cln_chemistry_.*\\.rda$", full.names = TRUE)
  clin_chemistry = lapply(chem_files, function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], e) })

  # These are reported in absolute terms and are written out that way; no log2 transform.
  # generate_clinical_DA.R applies its own log2 before fitting, so the modelling is unaffected.
  # There is one missing value; we treat it as missing completely at random.
  clin_chem_df = clin_chemistry %>%
    dplyr::bind_rows() %>%
    tibble::column_to_rownames("analyte_name")

  cln_out = clin_chem_df %>% tibble::rownames_to_column("feature_id") %>%
    dplyr::select(feature_id, dplyr::everything())

  # The analytes split across two omes by molecule class, so they are annotated and written
  # separately. .clinical_ome_for_feature() is the single source of truth for that routing.
  feature_omes = .clinical_ome_for_feature(cln_out$feature_id)

  prot_annotation = .annotate_clinical_prot(cln_out$feature_id[feature_omes == "prot-clinical"]) %>%
    dplyr::mutate(assay = "prot-clinical")
  # metab carries assay = "metab" with the ome in a platform column, matching the other metab
  # platforms (generate_metab_qc_norm.R writes assay = "metab" for every one of them); the prot
  # platforms instead name themselves in assay, so prot-clinical follows prot.
  # is_named / num_NAs / pct_na_imputed use the other metab platforms' column names and order.
  # Every analyte is named. Nothing is imputed here, so pct_na_imputed is a note rather than a
  # percentage.
  metab_platform_columns = cln_out %>%
    dplyr::transmute(feature_id,
                     is_named = TRUE,
                     num_NAs = rowSums(is.na(dplyr::across(-feature_id))),
                     pct_na_imputed = "No imputation was performed")
  metab_annotation = .annotate_clinical_metab(cln_out$feature_id[feature_omes == "metab-t-clinical"]) %>%
    dplyr::mutate(assay = "metab", platform = "metab-t-clinical") %>%
    dplyr::left_join(metab_platform_columns, by = "feature_id")
  all_annotation = dplyr::bind_rows(prot_annotation, metab_annotation)

  feature_metadata = cln_out %>%
    dplyr::select(feature_id) %>%
    dplyr::left_join(all_annotation, by = "feature_id")

  sample_metadata = pheno$data %>%
    dplyr::filter(vialLabel %in% colnames(cln_out))

  for(ome in unique(feature_omes)){
    ome_dir = file.path(freeze_output, .freeze_subdir_for_ome(ome))
    qc_output = file.path(ome_dir, "qc-norm")
    metadata_output = file.path(ome_dir, "metadata")
    dir.create(qc_output, recursive = TRUE, showWarnings = FALSE)
    dir.create(metadata_output, recursive = TRUE, showWarnings = FALSE)

    ome_features = cln_out$feature_id[feature_omes == ome]

    ome_qc_norm = cln_out %>%
      dplyr::filter(feature_id %in% ome_features)

    # the gene columns are all NA for the metabolites and the refmet columns are all NA for
    # the proteins; dropping them makes each file match its neighbours in the same directory.
    ome_feature_metadata = feature_metadata %>%
      dplyr::filter(feature_id %in% ome_features) %>%
      dplyr::select(assay, feature_id, dplyr::where(~!all(is.na(.x)))) %>%
      dplyr::relocate(dplyr::any_of("platform"), .after = assay)

    write_with_path_name(ome_qc_norm,
                         local_path = qc_output,
                         ome = ome,
                         tissue = "blood",
                         data_category = "qc-norm",
                         data_details = "absolute")

    write_with_path_name(ome_feature_metadata,
                         local_path = metadata_output,
                         ome = ome,
                         tissue = "blood",
                         data_category = "metadata",
                         data_details = "features")

    write_with_path_name(sample_metadata,
                         local_path = metadata_output,
                         ome = ome,
                         tissue = "blood",
                         data_category = "metadata",
                         data_details = "samples")
  }
}

.annotate_clinical_prot = function(prot_analyte_names){
  #analyte names don't match Ensembl external_gene_name; map manually to gene symbols first
  prot_analyte_to_gene = c("Insulin" = "INS", "Glucagon" = "GCG", "CK" = "CKM")

  biomart_anno = ensembl_v105_getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "uniprotswissprot"),
    filters = "external_gene_name",
    values = prot_analyte_to_gene[prot_analyte_names]
  ) %>%
    dplyr::distinct() %>%
    dplyr::filter(uniprotswissprot != "") %>%
    dplyr::mutate(ensembl_gene = str_remove(ensembl_gene_id, "\\..*"),
                  entrez_gene = as.character(entrezgene_id),
                  gene_symbol = external_gene_name,
                  feature_id = names(prot_analyte_to_gene)[match(external_gene_name, prot_analyte_to_gene)]) %>%
    dplyr::select(feature_id, entrez_gene, gene_symbol, ensembl_gene, uniprot = uniprotswissprot)

  return(biomart_anno)
}

.annotate_clinical_metab = function(metab_analyte_names){
  mets = paste(metab_analyte_names, collapse = "\n")

  h = curl::new_handle()
  curl::handle_setform(h, metabolite_name = mets)
  req = curl::curl_fetch_memory(
    "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_minID.php",
    handle = h
  )

  refmet_result = utils::read.csv(
    text = rawToChar(req$content),
    header = TRUE,
    na.strings = "-",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    sep = "\t"
  )

  refmet_result[is.na(refmet_result)] = '-'
  refmet_result[refmet_result == ''] = '-'

  metab_df = data.frame(feature_id = metab_analyte_names,
                        lookup_refmet = metab_analyte_names)

  refmet_output = refmet_result %>%
    dplyr::transmute(
      lookup_refmet = Input.name,
      refmet_name_std = Standardized.name,
      refmet_id = RefMet_ID,
      kegg_id = KEGG_ID
    ) %>%
    dplyr::left_join(metab_df, by = "lookup_refmet", relationship = "many-to-many") %>%
    dplyr::transmute(
      feature_id,
      refmet_name = dplyr::na_if(refmet_name_std, "-"),
      refmet_id = dplyr::na_if(refmet_id, "-"),
      kegg_id = dplyr::na_if(kegg_id, "-")
    ) %>%
    dplyr::distinct()

  return(refmet_output)
}

out_base <- .STAGING
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
generate_clinical_qc_norm(out_base)
message("Clinical chemistry freeze files written under ", file.path(out_base, "freeze"))
