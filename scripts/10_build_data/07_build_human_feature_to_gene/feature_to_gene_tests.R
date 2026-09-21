#!/usr/bin/env Rscript
# Tests for HUMAN_FEATURE_TO_GENE (step 07). Downstream: load_differential_analysis
# merges DA on key c("assay","feature_id"); run_cameraPR maps feature_id->gene_symbol
# (and prot-ph feature_id->flanking_sequence). Contract: keyed data.table, factor cols.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
suppressWarnings(suppressMessages(library(data.table)))
cat("== feature_to_gene tests ==\n")

hfg <- load_built("HUMAN_FEATURE_TO_GENE")
assert("HUMAN_FEATURE_TO_GENE built", !is.null(hfg), "data/HUMAN_FEATURE_TO_GENE.rda missing")
if (!is.null(hfg)) {
  assert("is a data.table", is.data.table(hfg))
  assert("key is c(assay, feature_id)", identical(key(hfg), c("assay", "feature_id")),
         sprintf("key = [%s]", paste(key(hfg), collapse = ", ")))
  assert("first two columns are assay, feature_id",
         identical(colnames(hfg)[1:2], c("assay", "feature_id")))
  assert_cols("HUMAN_FEATURE_TO_GENE", hfg, c("assay", "feature_id", "gene_symbol"))
  assert("nrow > 1e6", nrow(hfg) > 1e6, sprintf("nrow = %d", nrow(hfg)))
  assert("assay is character", is.character(hfg$assay))
  assert("gene_symbol is factor", is.factor(hfg$gene_symbol))
  assert("no duplicate (assay,feature_id) keys",
         !any(duplicated(hfg[, .(assay, feature_id)])))
  valid_assays <- c("epigen-atac-seq", "epigen-methylcap-seq", "metab", "metab-t-clinical",
                    "prot-ol", "prot-ph", "prot-pr", "prot-clinical", "transcript-rna-seq")
  assert_subset("assay values valid", unique(as.character(hfg$assay)), valid_assays)
  # run_cameraPR (step 12) keys phosphosite enrichment on the flanking sequence, not the
  # feature_id — PhosphoSitePlus kinase sets are defined over flanking sequences. The step-07
  # select used to drop this column, which silently reduced step 12 to the non-phospho omes;
  # it is now retained and asserted rather than skipped.
  assert("flanking_sequence present", "flanking_sequence" %in% colnames(hfg),
         "dropped by the step-07 select — prot-ph enrichment in step 12 needs it")
  if ("flanking_sequence" %in% colnames(hfg)) {
    ph <- as.character(hfg$flanking_sequence[as.character(hfg$assay) == "prot-ph"])
    assert("flanking_sequence populated for prot-ph", length(ph) > 0 && mean(!is.na(ph)) > 0.5,
           sprintf("%.1f%% non-NA over %d prot-ph row(s)", 100 * mean(!is.na(ph)), length(ph)))
  }

  # confident_site is deliberately absent. It is measured PER TISSUE and this object has no
  # tissue column, so carrying it would mean either a cross-tissue collapse (a value that is
  # not the measurement) or duplicate (assay, feature_id) keys for the 859 sites muscle and
  # adipose disagree on, which makes the map a one-to-many join for step 12. The exact
  # per-tissue value lives on *_PROT_PH_QC$feature_metadata, which is what step 15 filters
  # PTM-SEA input on.
  assert("confident_site absent", !"confident_site" %in% colnames(hfg),
         "per-tissue flag; read *_PROT_PH_QC$feature_metadata instead")

  # custom_annotation / relationship_to_gene -- the two peak-annotation columns the step-06
  # epigen stems write onto the metadata_features. The step-07 select used to drop them, so
  # the map carried an ATAC/methylcap feature to a gene without saying what the relationship
  # to that gene is: a promoter peak and one 40kb into an intron were indistinguishable.
  # Both are documented columns of this object (carry/04_document.R), so they are asserted
  # rather than skipped.
  epigen <- c("epigen-atac-seq", "epigen-methylcap-seq")
  is_epi <- as.character(hfg$assay) %in% epigen
  assert("custom_annotation present", "custom_annotation" %in% colnames(hfg),
         "dropped by the step-07 select -- epigen features lose their region class")
  assert("relationship_to_gene present", "relationship_to_gene" %in% colnames(hfg),
         "dropped by the step-07 select -- epigen features lose their distance to gene")

  if ("custom_annotation" %in% colnames(hfg)) {
    ca <- hfg$custom_annotation
    assert("custom_annotation is factor", is.factor(ca), sprintf("got %s", class(ca)[1]))
    assert("custom_annotation populated for epigen",
           sum(is_epi) > 0 && !anyNA(ca[is_epi]),
           sprintf("%d NA over %d epigen row(s)", sum(is.na(ca[is_epi])), sum(is_epi)))
    # Scoped to the epigen assays: a non-NA value elsewhere would mean bind_rows matched a
    # peak annotation onto a feature it does not describe.
    assert("custom_annotation NA outside epigen", all(is.na(ca[!is_epi])),
           sprintf("%d non-NA on non-epigen row(s)", sum(!is.na(ca[!is_epi]))))
    # ChIPseeker region classes, renamed by pre_cawg_get_peak_annotations_hs() in
    # lib/qc_helpers.R. A value outside this set means that renaming drifted.
    assert_subset("custom_annotation values valid",
                  unique(as.character(ca[is_epi])),
                  c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)",
                    "5' UTR", "3' UTR", "Exon", "Intron", "Overlaps Gene",
                    "Upstream (<5kb)", "Downstream (<5kb)", "Distal Intergenic"))
  }

  if ("relationship_to_gene" %in% colnames(hfg)) {
    rtg <- hfg$relationship_to_gene
    # Numeric, NOT the character it is read as or the factor the blanket as.factor() cast
    # would produce: the column exists to be compared (abs(x) < 5000, x > 0), and as a string
    # "40000" < "5000" is TRUE.
    assert("relationship_to_gene is numeric", is.numeric(rtg), sprintf("got %s", class(rtg)[1]))
    assert("relationship_to_gene populated for epigen",
           sum(is_epi) > 0 && !anyNA(rtg[is_epi]),
           sprintf("%d NA over %d epigen row(s)", sum(is.na(rtg[is_epi])), sum(is_epi)))
    assert("relationship_to_gene NA outside epigen", all(is.na(rtg[!is_epi])),
           sprintf("%d non-NA on non-epigen row(s)", sum(!is.na(rtg[!is_epi]))))
  }

  # Regression guard for the key. Both columns are coordinate-derived (computed from the peak
  # interval in the feature_id, not measured per tissue), so they need no cross-tissue
  # collapse -- but only as long as that holds. If it stops holding, the shared
  # feature_ids survive unique() as duplicate keys and the map becomes a one-to-many join.
  # ("no duplicate keys" above covers the whole object; this pins the cause to the epigen rows.)
  if (sum(is_epi) > 0) {
    for (a in epigen) {
      ia <- as.character(hfg$assay) == a
      if (sum(ia) > 0)
        assert_unique(sprintf("epigen annotation did not duplicate any %s key", a),
                      as.character(hfg$feature_id[ia]))
    }
  }
}
# cross-object: METABOLOMICS_CVS (step 05, upstream) re-derives refmet_name by joining
# HUMAN_FEATURE_TO_GENE metab rows, so every non-NA CV refmet must resolve here.
mc <- load_built("METABOLOMICS_CVS")
if (!is.null(hfg) && !is.null(mc)) {
  hmetab <- unique(as.character(hfg$refmet_name[as.character(hfg$assay) == "metab"]))
  rn <- unique(as.character(mc$refmet_name)); rn <- rn[!is.na(rn)]
  extra <- setdiff(rn, hmetab)
  # REPORT only: a mismatch here means the staged METABOLOMICS_CVS (step 05) is out of
  # sync with the rebuilt HFG (RefMet drift) — flag it, don't stop the build.
  report("METABOLOMICS_CVS refmet_name vs HUMAN_FEATURE_TO_GENE[metab]",
         if (length(extra) == 0) "all resolve"
         else sprintf("%d staged refmet not in rebuilt HFG, e.g. %s", length(extra),
                      paste(utils::head(extra, 3), collapse = ", ")))
} else skip("METABOLOMICS_CVS refmet_name vs HFG[metab]", "METABOLOMICS_CVS or HFG not built")

diff_vs_package("HUMAN_FEATURE_TO_GENE", hfg, tolerant = TRUE)
finish()
