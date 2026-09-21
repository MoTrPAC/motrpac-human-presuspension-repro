#!/usr/bin/env Rscript
# Stage 1 build step 14: UTORONTO_TFs  (deps: HUMAN_FEATURE_TO_GENE)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/UTORONTO_TFs.R.
#
# The phosphosites whose gene is a curated human transcription factor: the Human TFs
# database extract (Lambert et al. 2018) restricted to `Is TF? == "Yes"`, mapped through
# HUMAN_FEATURE_TO_GENE to every prot-ph feature sharing the gene symbol. Two columns,
# feature_id and gene_symbol.
#
# This object is a REGULATOR POOL, not an annotation table. `.load_scion_matrixes()` reads
# `UTORONTO_TFs$feature_id` to apply `subset_TFs = TRUE`, and `subset_qc()` guards its
# feature filter with `if (!is.null(desired_features))` — so an object without a
# `feature_id` column does not error, it silently turns the TF subsetting into a no-op and
# hands SCION every phosphosite in the matrix. That is the shape the released .rda is in
# (the raw 2,765 x 28 download), which is why the column set below is the whole contract.
#
# Output: scripts/10_build_data/data/UTORONTO_TFs.rda
#
# ---- Adaptations from upstream ----------------------------------------------------------
#
# 1. Reads the vendored sources/DatabaseExtract_v_1.01.txt instead of `~/Downloads/`, so the
#    build is offline and reproducible. See sources/README.md for provenance.
#
# 2. `assay`, not `platform`. HUMAN_FEATURE_TO_GENE has no `platform` column — its columns
#    are assay/feature_id/entrez_gene/gene_symbol/ensembl_gene/uniprot/refmet_name/
#    refmet_id/kegg_id/flanking_sequence. Upstream's `filter(platform == "prot-ph")` errors
#    against it, which is why that script cannot regenerate its own object as written.
#
# 3. The prot-ph filter runs BEFORE `distinct(feature_id)`, and `Is.TF.` before the join.
#    Upstream dedups first and filters after. A feature that maps to two genes gets a row
#    per gene, so deduping first keeps whichever row sorts first and can drop a real
#    prot-ph TF site whose other row is a non-TF gene or a different assay. Filtering first
#    makes the dedup a within-prot-ph tiebreak, which is what it is meant to be.
#
# 4. save() into scripts/10_build_data/data/ instead of usethis::use_data().

suppressWarnings(suppressMessages({ library(dplyr) }))

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT)
    ROOT <- dirname(ROOT)
}
HERE    <- file.path(ROOT, "scripts", "10_build_data", "14_build_utoronto_tfs")
out_dir <- file.path(ROOT, "scripts", "10_build_data", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a
load_one <- function(name) {
  f <- file.path(out_dir, paste0(name, ".rda"))
  if (!file.exists(f)) stop("missing ", name, ".rda — run its build step first: ", f)
  e <- new.env(); load(f, envir = e); e[[name]] %||% e[[ls(e)[1]]]
}

src <- file.path(HERE, "sources", "DatabaseExtract_v_1.01.txt")
if (!file.exists(src)) stop("missing the vendored TF database extract: ", src)

# check.names = TRUE is upstream's, and load-bearing: it is what makes the `Is TF?` and
# `HGNC symbol` headers addressable as Is.TF. and HGNC.symbol. [-1] drops the unnamed
# leading row-index column.
tfs <- read.csv(src, sep = "\t", check.names = TRUE)[-1] %>%
  dplyr::rename(gene_symbol = HGNC.symbol) %>%
  dplyr::filter(Is.TF. == "Yes") %>%
  dplyr::mutate(gene_symbol = as.character(gene_symbol))

HUMAN_FEATURE_TO_GENE <- load_one("HUMAN_FEATURE_TO_GENE")   # step 07
if (!"assay" %in% colnames(HUMAN_FEATURE_TO_GENE))
  stop("HUMAN_FEATURE_TO_GENE has no assay column — step 07 output is the wrong shape")

# Both keys forced to character: step 07 ships this object as a data.table whose
# gene_symbol and feature_id are factors over the full 1.9M-row level sets, and joining a
# factor to a character key coerces silently.
ftg <- HUMAN_FEATURE_TO_GENE %>%
  as.data.frame() %>%
  dplyr::filter(assay == "prot-ph") %>%
  dplyr::mutate(gene_symbol = as.character(gene_symbol),
                feature_id  = as.character(feature_id)) %>%
  dplyr::select(feature_id, gene_symbol)

UTORONTO_TFs <- tfs %>%
  dplyr::left_join(ftg, by = "gene_symbol", relationship = "many-to-many") %>%
  dplyr::filter(!is.na(feature_id)) %>%
  dplyr::distinct(feature_id, .keep_all = TRUE) %>%
  dplyr::select(feature_id, gene_symbol) %>%
  dplyr::arrange(feature_id)

if (nrow(UTORONTO_TFs) == 0)
  stop("UTORONTO_TFs is empty — no prot-ph feature matched a TF gene symbol")

save(UTORONTO_TFs, file = file.path(out_dir, "UTORONTO_TFs.rda"), compress = "xz")
message(sprintf("UTORONTO_TFs: %d prot-ph sites across %d TF genes (of %d curated TFs) -> %s",
                nrow(UTORONTO_TFs), dplyr::n_distinct(UTORONTO_TFs$gene_symbol),
                dplyr::n_distinct(tfs$gene_symbol), out_dir))
