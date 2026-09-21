#!/usr/bin/env Rscript
# Stage 1 data object: HUMAN_FEATURE_TO_GENE  (deps: OME_TISSUE_CODE, QC_NORM)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/HUMAN_FEATURE_TO_GENE.R.
#
# The feature->gene map is aggregated from every `metadata_features` file. Unlike
# the upstream script (which lists the staging bucket and dl_read_gcp's each file),
# this reads the LOCAL metadata_features written by the generate_normalized_expression
# stems into staging/freeze/*/metadata/ — so the map is built from THIS run's
# regenerated qc-norm outputs, no bucket access required. Instead of
# usethis::use_data() the object is save()d as a local .rda in
# scripts/10_build_data/data/.
#
# Deps:
#  - QC_NORM: the metadata_features are produced by the qc-norm stems; this reads
#    their local staging/freeze output, so a completed QC_NORM run is required.
#  - OME_TISSUE_CODE (preflight): the ome x tissue catalog used by the completeness
#    gate below to assert every feature-metadata-bearing combo actually has a file
#    before anything is combined, so a silently-missing assay fails loudly instead
#    of dropping out of the map.
#
# NOTE: this data/HUMAN_FEATURE_TO_GENE.rda is staging — it will eventually be
# organized as the HUMAN_FEATURE_TO_GENE data object in the Analysis package.
suppressWarnings(suppressMessages({ library(dplyr); library(data.table) }))

# Root from the PRECOVID_ROOT config env (config/pipeline.env); fall back to an
# upward search for config/pipeline.env so standalone runs work.
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
}
out_dir <- file.path(ROOT, "scripts", "10_build_data", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Freeze-side copy of the map: staging/freeze/resources/, mirroring where it sits on
# the bucket (resources/ is a top-level sibling of the per-ome-group folders, not
# nested under one of them).
resources_dir <- file.path(ROOT, "staging", "freeze", "resources")
dir.create(resources_dir, showWarnings = FALSE, recursive = TRUE)

# freeze_version(): per-file version from config/file_versions.json. Sourced directly
# rather than via lib/qc_helpers.R, which this step deliberately does not depend on
# (qc_helpers eagerly loads pheno + the other step-02/05 objects, none of which this
# stem needs). See lib/file_versions.R.
source(file.path(ROOT, "scripts", "10_build_data", "lib", "file_versions.R"))

# Upstream object: OME_TISSUE_CODE from preflight (the ome x tissue catalog).
otc_path <- file.path(ROOT, "scripts", "00_preflight", "data", "OME_TISSUE_CODE.rds")
if (!file.exists(otc_path))
  stop("missing OME_TISSUE_CODE.rds — build the 00_preflight data objects first: ", otc_path)
OME_TISSUE_CODE <- readRDS(otc_path)

# Local qc-norm output: the metadata_features written by generate_normalized_expression.
freeze_dir <- file.path(ROOT, "staging", "freeze")
feature_metadata_files <- list.files(freeze_dir, pattern = "metadata_features",
                                     recursive = TRUE, full.names = TRUE)
if (length(feature_metadata_files) == 0)
  stop("no metadata_features files under ", freeze_dir,
       " — run the QC_NORM stems (run_qc_norm.sh) first")

# --- Completeness gate ------------------------------------------------------
# Assert every feature-metadata-bearing ome x tissue combo has a file. Clinical
# chemistry (lab-*) carries no per-feature metadata, so it is excluded from the
# expectation.
expected <- OME_TISSUE_CODE %>%
  dplyr::filter(!grepl("^lab-", ome)) %>%
  dplyr::mutate(token = paste0(tissue_code, "_", ome))
have_file <- vapply(expected$token,
                    function(tk) any(grepl(tk, feature_metadata_files, fixed = TRUE)),
                    logical(1))
if (any(!have_file)) {
  missing <- expected[!have_file, c("tissue", "tissue_code", "ome")]
  stop("metadata_features file missing for ", sum(!have_file), " ome x tissue combo(s):\n",
       paste(sprintf("  %s / %s / %s", missing$tissue, missing$tissue_code, missing$ome),
             collapse = "\n"))
}
message("feature-metadata completeness: all ", nrow(expected),
        " expected ome x tissue combos have a metadata_features file")

# --- Aggregate into the feature -> gene map ---------------------------------
all_feat_to_gene <- list()
for (file in feature_metadata_files) {
  current_file <- read.csv(file, sep = "\t", colClasses = "character",
                           check.names = FALSE, stringsAsFactors = FALSE)
  all_feat_to_gene[[file]] <- current_file
}

HUMAN_FEATURE_TO_GENE <- dplyr::bind_rows(all_feat_to_gene) %>%
  dplyr::arrange(assay, feature_id) %>%
  dplyr::filter(is_named != "FALSE" | is.na(is_named)) %>%
  # any_of("flanking_sequence"): only the prot-ph metadata_features files carry it, so after
  # bind_rows it is NA for every other assay — and absent entirely if no prot-ph file is in
  # the freeze. It is kept because run_cameraPR() (step 12) keys phosphosite enrichment on
  # the flanking sequence rather than the feature_id: PhosphoSitePlus kinase sets (the PSP
  # database) are defined over flanking sequences. Dropping it here silently reduced step 12
  # to the four non-phospho omes.
  # The prot-ph and prot-pr per-tissue columns are deliberately NOT carried: confident_site,
  # confident_score, ptm_score and redundant_ids (prot-ph), num_peptides, percent_coverage,
  # protein_score and redundant_ids (prot-pr). This object is keyed on (assay, feature_id) with
  # no tissue column, so any value here would be a collapse across tissues rather than the
  # measurement -- muscle and adipose disagree on 15,825 of 20,664 shared sites for
  # confident_score, 5,495 of 6,079 shared proteins for num_peptides, and 2,162 / 539 for
  # redundant_ids. The select() below is a whitelist, so they are dropped by omission; the
  # duplicate-key guard would fail the build if one ever reached this point. Read them off
  # *_PROT_PH_QC / *_PROT_PR_QC $feature_metadata (step 08); confident_site is what step 15
  # filters PTM-SEA input on.
  # any_of("custom_annotation") / any_of("relationship_to_gene"): the two peak-annotation
  # columns .annotate_atac_features() and .annotate_methylcap_features() (step 06) write onto
  # the epigen metadata_features, so they are NA for every non-epigen assay after bind_rows and
  # absent entirely if no epigen file is in the freeze. They were dropped by this select even
  # though the object's own documentation (30_update_relevant_packages/carry/04_document.R)
  # describes both, so the map carried an ATAC/methylcap feature to a gene without saying what
  # the relationship to that gene IS -- a promoter peak and a peak 40kb into an intron mapped
  # to the same gene, indistinguishable. custom_annotation is the ChIPseeker region class
  # ("Promoter (<=1kb)", "Intron", "Distal Intergenic", ...) and relationship_to_gene the
  # signed base-pair distance to it, so a consumer can filter epigen features by proximity
  # without reaching back into the per-tissue QC objects.
  #
  # These need no cross-tissue collapse: both are derived from the peak COORDINATES (the
  # feature_id itself) via pre_cawg_get_peak_annotations_hs(), not measured per tissue, so
  # they agree on every feature_id shared between tissues -- 0 of the 306,788
  # ATAC and 0 of the 1,545,930 methylcap shared keys disagree in the current freeze. The
  # duplicate-key guard below pins that assumption rather than trusting it.
  dplyr::select(assay, feature_id, entrez_gene, gene_symbol, ensembl_gene,
                uniprot, refmet_name, refmet_id, kegg_id,
                dplyr::any_of("custom_annotation"),
                dplyr::any_of("relationship_to_gene"),
                dplyr::any_of("flanking_sequence"))

setDT(HUMAN_FEATURE_TO_GENE)

# --- relationship_to_gene: to numeric ---------------------------------------
# Read as character (colClasses = "character" above) and documented as numeric. It is already
# in the as.factor() exclusion list below, so without this cast it would silently ship as a
# character column and every comparison the column exists for -- abs(x) < 5000, x > 0 -- would
# either error or compare lexicographically ("40000" < "5000" is TRUE as a string).
if ("relationship_to_gene" %in% colnames(HUMAN_FEATURE_TO_GENE)) {
  rtg_chr <- trimws(as.character(HUMAN_FEATURE_TO_GENE[["relationship_to_gene"]]))
  rtg_chr[rtg_chr == ""] <- NA_character_
  rtg_num <- suppressWarnings(as.numeric(rtg_chr))
  bad <- !is.na(rtg_chr) & is.na(rtg_num)
  if (any(bad))
    stop("HUMAN_FEATURE_TO_GENE: relationship_to_gene has non-numeric value(s), e.g. ",
         paste(utils::head(unique(rtg_chr[bad]), 5), collapse = ", "))
  HUMAN_FEATURE_TO_GENE[, relationship_to_gene := rtg_num]
}

HUMAN_FEATURE_TO_GENE <- unique(HUMAN_FEATURE_TO_GENE)

# --- key uniqueness guard ---------------------------------------------------
# unique() drops rows identical across ALL columns, so a column that disagrees between two
# tissues for one feature_id survives as TWO rows with the same (assay, feature_id) key, which
# makes the map a one-to-many join for every consumer (step 12 joins DA results against it).
# The epigen annotation columns are asserted rather than collapsed (they are coordinate-derived
# and cannot legitimately differ), so if that ever stops holding the build must fail here
# instead of quietly double-counting every affected feature downstream.
dup_key <- duplicated(HUMAN_FEATURE_TO_GENE, by = c("assay", "feature_id"))
if (any(dup_key)) {
  ex <- HUMAN_FEATURE_TO_GENE[dup_key][seq_len(min(3L, sum(dup_key)))]
  stop("HUMAN_FEATURE_TO_GENE: ", sum(dup_key), " duplicate (assay, feature_id) key(s) -- ",
       "some column disagrees across tissues for the same feature and needs collapsing ",
       "or dropping. e.g.\n",
       paste(sprintf("  %s / %s", ex$assay, ex$feature_id), collapse = "\n"))
}

cols <- colnames(HUMAN_FEATURE_TO_GENE)
char_cols <- vapply(cols, class, character(1L)) == "character"
char_cols <- setdiff(names(char_cols), c("relationship_to_gene", "assay"))
HUMAN_FEATURE_TO_GENE[, (char_cols) := lapply(.SD, as.factor), .SDcols = char_cols]

setcolorder(x = HUMAN_FEATURE_TO_GENE, neworder = c("assay", "feature_id"))
setkeyv(x = HUMAN_FEATURE_TO_GENE, cols = c("assay", "feature_id"))

save(HUMAN_FEATURE_TO_GENE, file = file.path(out_dir, "HUMAN_FEATURE_TO_GENE.rda"),
     compress = TRUE, version = 3)
message("HUMAN_FEATURE_TO_GENE: ", nrow(HUMAN_FEATURE_TO_GENE), " rows -> ",
        file.path(out_dir, "HUMAN_FEATURE_TO_GENE.rda"))

# --- freeze copy as a resources .txt ----------------------------------------
# The bucket carries this map as a flat tab-separated table under resources/, so it
# is emitted here too rather than only as an .rda. Name and column set follow the
# staging bucket's v1.4 file, which is what this pipeline targets:
#   assay, feature_id, entrez_gene, gene_symbol, ensembl_gene, uniprot,
#   refmet_name, refmet_id, kegg_id
# plus the three annotation columns this pipeline keeps beyond that set:
# flanking_sequence (prot-ph only, step 12 keys phosphosite enrichment on it), and
# custom_annotation and relationship_to_gene (epigen only), each NA for every other
# assay.
# Written from the object, so they appear here automatically.

# Written UNQUOTED, matching the convention write_with_path_name() uses for every
# other freeze file — and with the same writer, write.table.
ftg_stem <- "motrpac-mappings-human-feature-to-gene"
ftg_txt  <- file.path(resources_dir,
                      paste0(ftg_stem, "_v", freeze_version(ftg_stem), ".txt"))
utils::write.table(as.data.frame(HUMAN_FEATURE_TO_GENE), ftg_txt, sep = "\t",
                   quote = FALSE, na = "NA", row.names = FALSE)
message("HUMAN_FEATURE_TO_GENE: ", nrow(HUMAN_FEATURE_TO_GENE), " rows -> ", ftg_txt)
