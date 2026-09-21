#!/usr/bin/env Rscript
# One-time builder for the local Ensembl v105 annotation cache.
#
# This is the ONLY script in Stage 1 that contacts Ensembl. It builds two
# artifacts into this folder so every QC-norm stem can annotate features OFFLINE
# and deterministically (the live Ensembl connection is flaky — intermittent
# curl_fetch_memory / "Fetch exons and CDS ... Error: [0]" / mirror-site errors —
# and must not be hit at pipeline run time):
#
#   1. txdb_hsapiens_ensembl_v105.sqlite   — TxDb for ChIPseeker peak annotation
#        (atac, methylcap), via txdbmaker::makeTxDbFromEnsembl + AnnotationDbi::saveDb.
#   2. ensembl_v105_gene_attributes.rds    — the full v105 gene-attribute table
#        (ensembl_gene_id, entrezgene_id, external_gene_name, uniprotswissprot),
#        the UNION of every getBM() attribute set across the stems. All 7 getBM
#        call sites are served by filtering this one table.
#
# Ensembl is flaky, so both fetches are retried a bounded number of times; if they
# still fail the script errors loudly (re-run when Ensembl is reachable).
#
# Each artifact is skipped if it already exists, so re-running is cheap and safe —
# `make sources` calls this on every run. Set FORCE=1 to rebuild both.
# See README.md for provenance.

RELEASE <- 105
OUT_DIR <- dirname(normalizePath(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
if (is.na(OUT_DIR) || !nzchar(OUT_DIR)) OUT_DIR <- getwd()
TXDB_PATH <- file.path(OUT_DIR, "txdb_hsapiens_ensembl_v105.sqlite")
GENE_PATH <- file.path(OUT_DIR, "ensembl_v105_gene_attributes.rds")

FORCE <- !Sys.getenv("FORCE") %in% c("", "0", "false", "FALSE")

MAX_TRIES <- 8
with_retry <- function(what, expr) {
  for (i in seq_len(MAX_TRIES)) {
    res <- tryCatch(force(expr), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    message(sprintf("[%s] attempt %d/%d failed: %s", what, i, MAX_TRIES, conditionMessage(res)))
    Sys.sleep(5)
  }
  stop(sprintf("%s: exhausted %d attempts (Ensembl unreachable?)", what, MAX_TRIES))
}

# --- 1. TxDb -----------------------------------------------------------------
if (file.exists(TXDB_PATH) && !FORCE) {
  message("TxDb cache already present, skipping build: ", TXDB_PATH)
} else {
  message("Building TxDb from Ensembl release ", RELEASE, " ...")
  txdb <- with_retry("makeTxDbFromEnsembl",
                     txdbmaker::makeTxDbFromEnsembl(organism = "Homo Sapiens", release = RELEASE))
  AnnotationDbi::saveDb(txdb, TXDB_PATH)
  message("  saved ", TXDB_PATH)
}

# --- 2. gene-attribute table -------------------------------------------------
# Requesting entrezgene_id + uniprotswissprot in ONE getBM makes BioMart drop
# every gene that has NEITHER xref (~44k non-coding/pseudogenes) — but the pinned
# stems query only 3 attributes and legitimately need those genes. So fetch each
# cross-reference SEPARATELY (2-attribute queries, keyed by ensembl_gene_id) and
# merge — this keeps all genes and reproduces each stem's narrower getBM exactly.
# Also batched (a no-filter getBM of all genes is too heavy for the v105 archive).
# This artifact is committed to git, so it is skipped unless absent or FORCE=1 —
# an unconditional refetch would rewrite the checked-in copy on every run.
if (file.exists(GENE_PATH) && !FORCE) {
  message("Gene-attribute table already present, skipping fetch: ", GENE_PATH)
} else {
  BATCH <- 2000
  message("Fetching v", RELEASE, " gene attributes (per-attribute, batched) ...")
  mart <- with_retry("useEnsembl",
                     biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = RELEASE))
  all_ids <- unique(with_retry("getBM gene-ids",
                               biomaRt::getBM(attributes = "ensembl_gene_id", mart = mart))$ensembl_gene_id)
  batches <- split(all_ids, ceiling(seq_along(all_ids) / BATCH))
  message("  ", length(all_ids), " genes in ", length(batches), " batches of ", BATCH)

  fetch_attr <- function(attr) {
    message("  fetching ", attr, " ...")
    parts <- lapply(seq_along(batches), function(i)
      with_retry(sprintf("getBM %s %d/%d", attr, i, length(batches)),
                 biomaRt::getBM(attributes = c("ensembl_gene_id", attr),
                                filters = "ensembl_gene_id", values = batches[[i]], mart = mart)))
    unique(do.call(rbind, parts))
  }
  sym <- fetch_attr("external_gene_name")
  ent <- fetch_attr("entrezgene_id")
  uni <- fetch_attr("uniprotswissprot")
  gene_tab <- Reduce(function(a, b) merge(a, b, by = "ensembl_gene_id", all = TRUE),
                     list(sym, ent, uni))
  gene_tab <- unique(gene_tab[, c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "uniprotswissprot")])
  saveRDS(gene_tab, GENE_PATH)
  message(sprintf("  saved %s (%d rows, %d unique genes)",
                  GENE_PATH, nrow(gene_tab), length(unique(gene_tab$ensembl_gene_id))))
}
message("Ensembl v", RELEASE, " cache build complete.")
