#!/usr/bin/env Rscript
# Stage 1 QC-norm stem: epigen-atac-seq (ATAC-seq)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_normalized_expression/
# generate_atac_qc_norm.R.
#
# In-house objects/functions (load_pheno, process_covariates, write_with_path_name)
# come from the vendored lib/qc_helpers.R. Raw qa-qc-metrics + peak-count paths are
# resolved from the preflight quant-id catalog via quantid_path() (upstream
# hard-coded human-precovid/results paths are dead). Tissue codes are mapped
# explicitly (muscle -> t06-muscle; the t10-muscle-powder ATAC is a reference
# standard and is excluded; blood -> t05-pbmc). Gated dl_read_gcp + variancePartition
# voom/dream + limma, and the biomaRt/ChIPseeker/txdbmaker peak annotation (online),
# are unchanged. Adds the metadata<->matrix alignment before normalization, and drops the
# OUTLIERS samples from the counts the way every other qc-norm stem does.
#
# RERUN_ATAC (config/pipeline.env) selects between two paths. TRUE, the default, is the
# adapted upstream build described above. FALSE calls .stage_atac_qc_norm() instead, which
# copies the released qc-norm matrix into the freeze verbatim and regenerates only the
# metadata — the same treatment generate_methylcap_qc_norm.R gives MethylCap-seq.
#
# Output: BIC-named freeze files under staging/freeze/epigenomics/.
suppressWarnings(suppressMessages({ library(dplyr); library(tibble); library(stringr); library(data.table) }))
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

# epigen-atac-seq tissue -> quant-id tissue_code (t10-muscle-powder excluded: reference standard)
.ATAC_TC <- c(muscle = "t06-muscle", blood = "t05-pbmc")

# PARALLEL PATH: voomWithDreamWeights is variancePartition, so it shares the step-09 dream
# fits' backend — now MulticoreParam (fork), which does not lose features. Under the SOCK
# backend it replaced, one feature per worker vanished and the fit errored outright when
# features <= workers: 100 peaks -> 90 on 10 workers; 20 -> 10; 10, 5 and 3 all error. Default
# comes from variancepartition_parallel_enabled() (VARIANCEPARTITION_PARALLEL_CORES in
# config/pipeline.env), not from PARALLEL_CORES, which still drives mice.
generate_atac_qc_norm = function(repo_local_dir, parallel = variancepartition_parallel_enabled(),
                                 rerun = rerun_atac_enabled()){
  # RERUN_ATAC=FALSE stages the released matrix instead of fitting one. Branch before
  # the variancePartition check: that path needs none of the modelling stack.
  if (!rerun) return(invisible(.stage_atac_qc_norm(repo_local_dir)))
  if (!requireNamespace("variancePartition", quietly = TRUE))
    stop("variancePartition required for ATAC voom/dream normalization")
  desired_ome = 'epigen-atac-seq'; tissue_types = c('muscle', 'blood')
  local_path = repo_local_dir
  ome_meta = MotrpacBicQC::dl_read_gcp(quantid_path(desired_ome, data_details = "qa-qc-metrics"),
                                       tmpdir = .RAW_FILES, sep = ",")
  ome_meta$vialLabel = as.character(ome_meta$vialLabel)

  # Full pheno, as every other qc-norm stem uses: the ADU_PAS training samples (5 muscle,
  # 6 pbmc) stay in the matrix. The DA stem fits acute contrasts only and filters to
  # ADU_BAS itself (generate_atac_da.R).
  pheno_data_parsed = load_pheno(load_acute_only = FALSE)$pheno_data
  metadata_path = paste0(local_path, "freeze/epigenomics/metadata/")
  qc_norm_path = paste0(local_path, "freeze/epigenomics/qc-norm/")
  dir.create(metadata_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_norm_path, recursive = TRUE, showWarnings = FALSE)

  for (tissue in tissue_types){
    message(paste("Generating normalized matrixes for", tissue))
    file_load = quantid_path(desired_ome, tissue_code = .ATAC_TC[[tissue]], data_details = "counts")
    raw_atac_input = MotrpacBicQC::dl_read_gcp(file_load, sep = '\t', tmpdir = .RAW_FILES)

    raw_atac_input = raw_atac_input %>%
      dplyr::mutate(rownames = stringr::str_c(chrom,":",start,"-",end,sep = "")) %>%
      tibble::column_to_rownames("rownames") %>%
      dplyr::filter(!(chrom %in% c("chrX","chrY"))) %>%
      dplyr::select(-chrom,-start,-end) %>%
      as.data.frame()

    outliers_vialLabels = OUTLIERS$vialLabel
    desired_samples = setdiff(colnames(raw_atac_input), outliers_vialLabels) #counts minus outliers
    tissue_pheno = pheno_data_parsed[pheno_data_parsed$vialLabel %in% desired_samples, ]
    tissue_metadata = ome_meta[ome_meta$vialLabel %in% tissue_pheno$vialLabel, ]
    raw_atac_input = raw_atac_input[, colnames(raw_atac_input) %in% tissue_pheno$vialLabel] #remove HA, peds

    min_count = 2 * stats::median(as.matrix(raw_atac_input)) #aggressive pruning for ATAC
    min_samples = 0.5*dim(raw_atac_input)[2]

    atac_filtered = raw_atac_input[rowSums(data.frame(lapply(raw_atac_input, function(x) as.numeric(x >= min_count)), check.names=FALSE)) >= min_samples,]
    meta = merge(tissue_pheno, tissue_metadata, by = "vialLabel")
    # align metadata rows to the count-matrix columns (samples) via match() before
    # normalization/batch correction — merge() sorts by vialLabel, otherwise the
    # covariates would be applied to the wrong samples.
    atac_filtered = atac_filtered[, colnames(atac_filtered) %in% meta$vialLabel, drop = FALSE]
    meta = meta[match(colnames(atac_filtered), meta$vialLabel), ]
    meta = tibble::column_to_rownames(meta, "vialLabel")
    process_metadata = process_covariates(meta = meta,
                                          selected_ome = desired_ome,
                                          tissue_input = tissue)
    meta = process_metadata$metadata

    dge_list <- edgeR::DGEList(counts = atac_filtered)
    dge_list <- edgeR::calcNormFactors(dge_list)
    formula = process_metadata$full_formula
    message(paste(tissue, "atac", formula))

    if (parallel){
      # cores from config/pipeline.env (VARIANCEPARTITION_PARALLEL_CORES); fall back to detectCores()-2
      num_cores = suppressWarnings(as.integer(Sys.getenv("VARIANCEPARTITION_PARALLEL_CORES", "")))
      if (is.na(num_cores) || num_cores < 1) num_cores = parallel::detectCores() - 2
      param <- BiocParallel::MulticoreParam(num_cores, progressbar = TRUE)
      suppressWarnings({voom_object <- variancePartition::voomWithDreamWeights(dge_list, formula = stats::as.formula(formula), data = meta, BPPARAM = param)})
    }else{
      voom_object <- variancePartition::voomWithDreamWeights(dge_list, formula = stats::as.formula(formula), data = meta)
    }

    # Feature-count guard. Kept after the move to MulticoreParam: the SOCK loss it was written
    # for is gone, but nothing downstream notices a short matrix, so the fork backend is
    # checked rather than trusted too. See VARIANCEPARTITION_PARALLEL_CORES in config/pipeline.env.
    if (nrow(voom_object$E) != nrow(dge_list))
      stop(sprintf("atac %s: voomWithDreamWeights returned %d of %d features (parallel = %s) — %d lost",
                   tissue, nrow(voom_object$E), nrow(dge_list), parallel,
                   nrow(dge_list) - nrow(voom_object$E)))

    atac_norm <- voom_object$E

    technical_cov = paste(process_metadata[["technical_cov"]]$covariate, collapse = " + ")
    design_cov = paste(process_metadata[["design_cov"]], collapse = " + ")
    message(tissue," technical: ", technical_cov, " design: ", design_cov)
    batch_corrected = limma::removeBatchEffect(atac_norm,
                                               covariates = stats::model.matrix(stats::as.formula(paste("~ ", technical_cov)), data = meta),
                                               design = stats::model.matrix(stats::as.formula(paste("~ ", design_cov)), data = meta))
    batch_corrected = as.data.frame(batch_corrected)
    batch_corrected$feature_id = rownames(batch_corrected)
    batch_corrected = batch_corrected %>% dplyr::select(feature_id, everything())

    only_features = batch_corrected %>% dplyr::select(feature_id)
    atac_annotated = .annotate_atac_features(only_features)

    write_with_path_name(tissue_metadata, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'samples')
    write_with_path_name(batch_corrected, local_path = qc_norm_path, ome = desired_ome, tissue = tissue, data_category = 'qc-norm', data_details = 'log-cpm')
    #in version 1.4 we also attach the gene level information into the feature metadata.
    write_with_path_name(atac_annotated, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'features')
  }
}

# RERUN_ATAC=FALSE path. Does for ATAC what generate_methylcap_qc_norm.R does for
# methylcap: the released qc-norm matrix is treated as a vendored INPUT and copied into
# the freeze verbatim, and is consumed only for its feature_ids and sample columns. The
# feature metadata (gene annotation) and sample metadata are still generated here, by the
# same code and from the same sources as the rerun path.
#
# The vendored filename already carries the standard BIC
# <tissue_code>_<ome>_qc-norm_log-cpm_v<x.y> token, so once copied it sits alongside every
# other ome's qc-norm and the freeze tests, step 08 and everything downstream pick it up
# with no special case.
.stage_atac_qc_norm = function(repo_local_dir){
  desired_ome = 'epigen-atac-seq'
  local_path = repo_local_dir
  metadata_path = paste0(local_path, "freeze/epigenomics/metadata/")
  qc_norm_path = paste0(local_path, "freeze/epigenomics/qc-norm/")
  dir.create(metadata_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_norm_path, recursive = TRUE, showWarnings = FALSE)

  ome_meta = MotrpacBicQC::dl_read_gcp(quantid_path(desired_ome, data_details = "qa-qc-metrics"),
                                       tmpdir = .RAW_FILES, sep = ",")
  ome_meta$vialLabel = as.character(ome_meta$vialLabel)
  pheno_data_parsed = load_pheno(load_acute_only = FALSE)$pheno_data

  # See sources/atac_qc_norm/../README.md for provenance. Presence is enforced by
  # scripts/10_build_data/check_required_inputs.sh when RERUN_ATAC=FALSE.
  src_dir = file.path(.QC_REPO, "scripts", "10_build_data", "06_generate_qc_norm",
                      "sources", "atac_qc_norm")
  src_files = list.files(src_dir, pattern = "atac.*qc-norm.*\\.txt$", full.names = TRUE)

  # Driven by the expected tissue list rather than by whatever is in src_dir, so a
  # missing or duplicated tissue is an error instead of a silently shorter freeze.
  for (tissue in names(.ATAC_TC)){
    tissue_code = .ATAC_TC[[tissue]]
    hit = src_files[grepl(tissue_code, basename(src_files), fixed = TRUE)]
    if (length(hit) != 1)
      stop("expected exactly 1 vendored ATAC qc-norm file for ", tissue_code,
           " but found ", length(hit), " in ", src_dir)
    message(paste("Staging released qc-norm for", tissue))

    dest = file.path(qc_norm_path, basename(hit))
    if (!file.copy(hit, dest, overwrite = TRUE))
      stop("copy failed: ", hit, " -> ", dest)
    if (file.size(dest) != file.size(hit))
      stop("size mismatch after copying ", basename(hit), " — source ", file.size(hit),
           " bytes, freeze copy ", file.size(dest), " bytes")

    qc_norm = read.csv(hit, sep = "\t", check.names = FALSE)
    # first column holds the feature_id (chr:start-end); the rest are vialLabels
    only_features = data.frame(feature_id = qc_norm[[1]], stringsAsFactors = FALSE)
    atac_annotated = .annotate_atac_features(only_features)

    # sed-adult cohort scoping, same as the rerun path: pheno defines who is in scope
    # and the qa-qc table is then subset to those vials.
    sample_cols = setdiff(colnames(qc_norm), colnames(qc_norm)[1])
    tissue_pheno = pheno_data_parsed[pheno_data_parsed$vialLabel %in% sample_cols, ]
    tissue_metadata = ome_meta[ome_meta$vialLabel %in% tissue_pheno$vialLabel, ]
    message(sprintf("  %s: %d released sample(s), %d in pheno, %d with qa-qc metrics",
                    tissue, length(sample_cols), nrow(tissue_pheno), nrow(tissue_metadata)))

    write_with_path_name(tissue_metadata, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'samples')
    write_with_path_name(atac_annotated, local_path = metadata_path, ome = desired_ome, tissue = tissue, data_category = 'metadata', data_details = 'features')
  }
}

.annotate_atac_features = function(feature_metadata){
  atacpeakmeta = feature_metadata %>%
    dplyr::mutate(chrom = gsub(":.*", "", feature_id),
                  start = as.numeric(gsub(".*:|-.*", "", feature_id)),
                  end = as.numeric(gsub(".*-", "", feature_id))) %>%
    data.table::as.data.table()

  atac_peakdf = pre_cawg_get_peak_annotations_hs(atacpeakmeta)
  attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name")
  atac_lookup_df = ensembl_v105_getBM(attributes = attributes,
                                  filters = "ensembl_gene_id",
                                  values = atac_peakdf$ensembl_gene) %>%
    dplyr::full_join(atac_peakdf, c("ensembl_gene_id" = "ensembl_gene")) %>%
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
                  assay = "epigen-atac-seq") %>%
    dplyr::relocate(assay)

  return(atac_lookup_df)
}

# pre_cawg_get_peak_annotations_hs() now lives in lib/qc_helpers.R (shared with
# the methylcap stem); it is provided by the source() at the top of this file.

out_base <- .STAGING
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
generate_atac_qc_norm(paste0(out_base, "/"))
message(if (rerun_atac_enabled()) "ATAC QC-norm written under " else "ATAC QC-norm staged under ",
        file.path(out_base, "freeze", "epigenomics"))
