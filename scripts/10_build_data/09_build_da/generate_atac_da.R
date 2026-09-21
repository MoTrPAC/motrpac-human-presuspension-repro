#!/usr/bin/env Rscript
# Step 09 DA stem: epigen-atac-seq (blood/t05-pbmc and muscle/t06-muscle).
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_differential_analysis/
# generate_atac_DA.R.
#
# Like transcriptomics, ATAC is fit on RAW COUNTS with voom weighting: the qc-norm matrix
# supplies the peak and sample selection, the model runs on the counts. Technical covariates
# are included (frip), matching upstream.
#
# Upstream calls load_qc(..., epigen = TRUE) to reach the epigen objects shipped in the
# installed Data package. Step 08 now builds {TISSUE}_EPIGEN_ATAC_SEQ_QC locally, so this
# stem reads it through load_qc_local() exactly like every other ome — no epigen special
# case is needed on either side.
#
# RERUN_ATAC (config/pipeline.env) selects between two paths. TRUE, the default, is the
# adapted upstream fit described above. FALSE calls .stage_atac_da() instead, which copies
# the released DA tables into the freeze verbatim — the same treatment
# generate_methylcap_da.R gives MethylCap-seq, and it needs neither the dream engine nor
# the step 05-08 objects.
#
# Output: staging/freeze/epigenomics/da/
#   human-precovid-sed-adu_{t05-pbmc,t06-muscle}_epigen-atac-seq_da_dream-acute_v<x.y>.txt

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))
# da_common.R pulls in the dream engine and the whole modelling chain, which is only
# needed to FIT the model. Under RERUN_ATAC=FALSE this stem is a file copy, so it is left
# unsourced and the freeze stops being gated on steps 05-08 — the same reasoning
# generate_methylcap_da.R documents for skipping it outright.
if (rerun_atac_enabled())
  source(file.path(ROOT, "scripts", "10_build_data", "09_build_da", "da_common.R"))
suppressWarnings(suppressMessages({library(tibble); library(stringr); library(data.table)}))

# ---- Upstream build, verbatim from .../generate_differential_analysis/generate_atac_DA.R ----
#
# .generate_atac_inputs = function(repo_local_dir,
#                                  model_type,
#                                  tissue,
#                                  parallel = FALSE){
#   check_package_installation("edgeR")
#   message("Note: Epigen ATAC Seq only has 5 training samples, all of which are control. No training analysis will be done.")
#   desired_ome = 'epigen-atac-seq'
#   local_path =  file.path(repo_local_dir, "data/tmp/")
#   data_path = .find_path_name(desired_ome = desired_ome, tissue = tissue, data_type = "epigen-atac-seq_counts")
#   if (length(data_path) > 0 & model_type == "acute") {
#     raw_counts_input = MotrpacBicQC::dl_read_gcp(data_path, sep = '\t', tmpdir = local_path) %>%
#       dplyr::mutate(rownames = stringr::str_c(chrom,":",start,"-",end,sep = "")) %>%
#       tibble::column_to_rownames("rownames") %>%
#       dplyr::filter(!(chrom %in% c("chrX","chrY"))) %>%
#       dplyr::select(-chrom,-start,-end) %>%
#       as.data.frame()
#
#     parsed_qc_norm = load_qc(selected_omes = desired_ome,
#                              selected_tissues = tissue,
#                              epigen = TRUE,
#                              load_acute_only = FALSE)
#     #-> so here we want to make sure the genes we chose and the participants we chose are the same as the ones that pass through expression logcpm cutoffs
#     metadata = parsed_qc_norm[[tissue]][[desired_ome]][['sample_metadata']] %>%
#       dplyr::filter(visitcode == 'ADU_BAS')  #filter just to the initial acute bout
#
#     rownames(metadata) = metadata$vialLabel
#     raw_counts_input = raw_counts_input %>%
#       dplyr::filter(rownames(.) %in% rownames(parsed_qc_norm[[tissue]][[desired_ome]][['qc_norm']]))  %>% #filter to select genes
#       dplyr::select(as.character(metadata$vialLabel)) #reorganize and filter to only desired participants
#
#     pre_dge <- edgeR::DGEList(counts = raw_counts_input)
#     pre_dge <- edgeR::calcNormFactors(pre_dge)
#
#     process_metadata = process_covariates(meta = metadata,
#                                           selected_ome = desired_ome,
#                                           tissue_input = tissue)
#
#     .run_models(repo_local_dir = repo_local_dir,
#                 model_type = model_type,
#                 expression_object = pre_dge,
#                 process_metadata = process_metadata,
#                 tissue = tissue,
#                 ome = desired_ome,
#                 voom = TRUE,
#                 parallel = parallel)
#   }
# }
#
# Adapted below: load_qc(epigen = TRUE) becomes .atac_leaf_from_freeze() (see the note at the
# top), .find_path_name()'s live gsutil listing becomes quantid_path() over the preflight
# catalog, and .run_models() becomes run_da_models(). Modelling is unchanged.


# PARALLEL PATH: `parallel` reaches run_dream()'s MulticoreParam fork backend. See
# VARIANCEPARTITION_PARALLEL_CORES in config/pipeline.env for why it is fork and not SOCK.
generate_atac_da = function(tissue, parallel = da_parallel_enabled()){
  desired_ome = "epigen-atac-seq"
  message("Note: Epigen ATAC Seq only has 5 training samples, all of which are control. No training analysis will be done.")
  tissue_code = .match_ome_tissue_code(desired_ome = desired_ome, input_tissue = tissue)

  data_path = quantid_path(desired_ome, tissue_code = tissue_code, data_details = "counts")
  if (length(data_path) == 0)
    stop("no ATAC counts path in the quant-id catalog for ", tissue, "/", desired_ome)

  raw_counts_input = MotrpacBicQC::dl_read_gcp(data_path, sep = "\t",
                                               tmpdir = .RAW_FILES, check_first = TRUE) %>%
    dplyr::mutate(rownames = stringr::str_c(chrom, ":", start, "-", end, sep = "")) %>%
    tibble::column_to_rownames("rownames") %>%
    dplyr::filter(!(chrom %in% c("chrX", "chrY"))) %>%
    dplyr::select(-chrom, -start, -end) %>%
    as.data.frame()

  parsed_qc_norm = load_qc_local(selected_omes = desired_ome,
                                 selected_tissues = tissue,
                                 load_acute_only = FALSE)
  leaf = parsed_qc_norm[[tissue]][[desired_ome]]
  if (is.null(leaf))
    stop("no local *_QC object for ", tissue, "/", desired_ome, " — run step 08 first")

  metadata = leaf[["sample_metadata"]] %>% dplyr::filter(visitcode == "ADU_BAS")
  rownames(metadata) = metadata$vialLabel

  raw_counts_input = raw_counts_input %>%
    dplyr::filter(rownames(.) %in% rownames(leaf[["qc_norm"]])) %>%
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

# RERUN_ATAC=FALSE path. Does for the ATAC DA tier what generate_methylcap_da.R does for
# methylcap: the released tables are treated as a vendored INPUT and copied into the
# freeze verbatim, not re-fit. Their filenames already carry the standard BIC
# <tissue_code>_<ome>_da_dream-acute_v<x.y> token, so once copied they sit alongside every
# other ome's DA output with no special case downstream.
.stage_atac_da = function(){
  desired_ome = "epigen-atac-seq"
  # Same path da_common.R's .da_freeze_dir() builds, but that file is not sourced on this
  # path, so the shared .freeze_subdir_for_ome() from qc_helpers.R is used directly.
  da_path = file.path(.STAGING, "freeze", .freeze_subdir_for_ome(desired_ome), "da")
  dir.create(da_path, recursive = TRUE, showWarnings = FALSE)

  # See sources/atac_da/../README.md for provenance. Presence is enforced by
  # scripts/10_build_data/check_required_inputs.sh when RERUN_ATAC=FALSE.
  src_dir = file.path(ROOT, "scripts", "10_build_data", "09_build_da", "sources", "atac_da")
  src_files = list.files(src_dir, pattern = "atac.*_da_.*\\.txt$", full.names = TRUE)

  # Driven by the expected tissue list rather than by whatever is in src_dir, so a missing
  # or duplicated tissue is an error instead of a silently shorter freeze.
  for (tissue in c("blood", "muscle")) {
    tissue_code = .match_ome_tissue_code(desired_ome = desired_ome, input_tissue = tissue)
    hit = src_files[grepl(tissue_code, basename(src_files), fixed = TRUE)]
    if (length(hit) != 1)
      stop("expected exactly 1 vendored ATAC DA file for ", tissue_code,
           " but found ", length(hit), " in ", src_dir)

    # basename(hit), NOT freeze_version(): this path deliberately diverges from the map.
    # file_versions.json carries 2.0 for the ATAC DA tokens, which is the version
    # generate_atac_da() stamps when it actually FITS the model (run_da_models passes
    # version = NULL, so write_with_path_name resolves the map). Here nothing is fit — the
    # vendored table is copied through byte-for-byte — so it keeps the v1.2 it was released
    # under. Renaming an unmodified external artifact to 2.0 would assert a content change
    # that did not happen. The map entry is therefore the RERUN_ATAC=TRUE version, and a
    # v1.2 file in the freeze under RERUN_ATAC=FALSE is correct, not stale.
    dest = file.path(da_path, basename(hit))
    if (!file.copy(hit, dest, overwrite = TRUE))
      stop("copy failed: ", hit, " -> ", dest)
    # A truncated copy would still parse as a valid table, just with features missing,
    # so the size is compared rather than trusting file.copy()'s return alone.
    if (file.size(dest) != file.size(hit))
      stop("size mismatch after copying ", basename(hit), " — source ", file.size(hit),
           " bytes, freeze copy ", file.size(dest), " bytes")

    message(sprintf("%s / %s: copied %s (%.1f MB) -> %s",
                    tissue, desired_ome, basename(hit), file.size(dest) / 2^20, da_path))
  }
}

if (rerun_atac_enabled()) {
  for (tis in c("blood", "muscle")) generate_atac_da(tissue = tis)
} else {
  .stage_atac_da()
}
