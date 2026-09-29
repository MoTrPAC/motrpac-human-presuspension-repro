# Stage 3 carry, step 4 — bring the documentation back in line with what the
# packages now ship, then re-roxygenise.
#
# R CMD check does not read `\describe{}` against the data, so a stale @format is
# invisible to it: an object can gain and lose columns through every check a
# package has and still be documented as its previous self. The edits here are
# therefore driven off the carried objects rather than off a hand-kept list, and
# each one asserts its anchor before touching the file — an edit that silently
# matches nothing is the failure mode worth guarding against.
#
# Usage:
#   Rscript 04_document.R --manifest <tsv> --routing <tsv> --out-root <dir> --out <dir>
#
# Writes:
#   <out>/document_report.tsv

.here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(.here, "lib", "carry_helpers.R"))

opt <- parse_opts(commandArgs(trailingOnly = TRUE),
                  c("manifest", "routing", "out-root", "out", "data-pkg", "analysis-pkg"))
for (need in c("manifest", "routing", "out-root", "out")) {
  if (is.null(opt[[need]])) stop("missing required option --", need)
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
rpt <- new_report()

manifest <- read_tsv(opt$manifest)
routing  <- read_tsv(opt$routing)

PKG <- c(data = "MotrpacHumanPreSuspensionData",
         analysis = "MotrpacHumanPreSuspensionAnalysis")
root <- function(dest) file.path(opt$`out-root`, PKG[[dest]])

# ---- Editing primitives ------------------------------------------------------
# Both assert before they write, and both are idempotent, so a re-run of the step
# is a no-op rather than a second insertion.

#' Insert lines after the first line matching `anchor` exactly
insert_after <- function(path, anchor, new_lines, check, already = NULL) {
  lines <- readLines(path, warn = FALSE)
  if (!is.null(already) && any(lines == already)) {
    record(rpt, "PASS", check, "already present — no change")
    return(invisible(FALSE))
  }
  i <- which(lines == anchor)
  if (!length(i)) {
    record(rpt, "FAIL", check, paste0("anchor not found in ", basename(path), ": ", anchor))
    return(invisible(FALSE))
  }
  i <- i[1]
  writeLines(append(lines, new_lines, after = i), path)
  record(rpt, "PASS", check, sprintf("inserted %d line(s) after '%s' in %s",
                                     length(new_lines), anchor, basename(path)))
  invisible(TRUE)
}

#' Replace the inclusive line range from `from` to `to` with `new_lines`
replace_block <- function(path, from, to, new_lines, check) {
  lines <- readLines(path, warn = FALSE)
  i <- which(lines == from)
  if (!length(i)) {
    record(rpt, "FAIL", check, paste0("start anchor not found in ", basename(path), ": ", from))
    return(invisible(FALSE))
  }
  i <- i[1]
  j <- which(lines == to & seq_along(lines) > i)
  if (!length(j)) {
    record(rpt, "FAIL", check, paste0("end anchor not found in ", basename(path), ": ", to))
    return(invisible(FALSE))
  }
  j <- j[1]
  writeLines(c(lines[seq_len(i - 1L)], new_lines, lines[(j + 1L):length(lines)]), path)
  record(rpt, "PASS", check, sprintf("replaced %d line(s) in %s", j - i + 1L, basename(path)))
  invisible(TRUE)
}

#' Append a `@rdname` stanza for a data object to its family's file
append_stanza <- function(path, rdname, object, check) {
  lines <- readLines(path, warn = FALSE)
  quoted <- paste0('"', object, '"')
  if (any(lines == quoted)) {
    record(rpt, "PASS", check, paste0(object, " already documented"))
    return(invisible(FALSE))
  }
  writeLines(c(lines, "", paste0("#' @rdname ", rdname), "#' @format NULL",
               "#' @usage NULL", quoted), path)
  record(rpt, "PASS", check, paste0("added ", object, " to ", rdname))
  invisible(TRUE)
}

# ---- The objects new to the packages this cycle ------------------------------
# v1.3 carried clinical chemistry as one assay; v2.0 splits it into a metabolomics
# and a proteomics assay. Each half needs a doc entry in its family block, or
# R CMD check reports an undocumented data object.
#
# The three {TISSUE}_METAB_SUM_STATS objects are the other half of a rename: step
# 11 stacks the research metabolomics platforms per tissue the way the *_DA
# objects are keyed, so the per-platform objects leave through ORPHAN_DROP and
# their documentation with them. One spec per tissue, anchored on the family
# block's tissue heading rather than on a per-platform name, which the withdrawal
# pass below deletes.

NEW_OBJECTS <- list(
  list(dest = "data",     file = "R/QC_RESULTS.R",  rdname = "QC_RESULTS",
       anchor = "#' BLOOD_PROT_OL_QC",
       objects = c("BLOOD_METAB_T_CLINICAL_QC", "BLOOD_PROT_CLINICAL_QC")),
  list(dest = "analysis", file = "R/DA_RESULTS.R",  rdname = "DA_RESULTS",
       anchor = "#' BLOOD_PROT_OL_DA",
       objects = c("BLOOD_METAB_T_CLINICAL_DA", "BLOOD_PROT_CLINICAL_DA")),
  list(dest = "analysis", file = "R/SUM_STATS.R",   rdname = "SUM_STATS_RESULTS",
       anchor = "#' BLOOD_PROT_OL_SUM_STATS",
       objects = c("BLOOD_METAB_T_CLINICAL_SUM_STATS", "BLOOD_PROT_CLINICAL_SUM_STATS")),
  list(dest = "analysis", file = "R/SUM_STATS.R",   rdname = "SUM_STATS_RESULTS",
       anchor = "#' ## Adipose", objects = "ADIPOSE_METAB_SUM_STATS"),
  list(dest = "analysis", file = "R/SUM_STATS.R",   rdname = "SUM_STATS_RESULTS",
       anchor = "#' ## Blood",   objects = "BLOOD_METAB_SUM_STATS"),
  list(dest = "analysis", file = "R/SUM_STATS.R",   rdname = "SUM_STATS_RESULTS",
       anchor = "#' ## Muscle",  objects = "MUSCLE_METAB_SUM_STATS")
)

for (spec in NEW_OBJECTS) {
  path <- file.path(root(spec$dest), spec$file)
  if (!file.exists(path)) {
    record(rpt, "FAIL", paste0("doc:", basename(spec$file)), "file not found in the test package")
    next
  }
  insert_after(path, spec$anchor, paste0("#' ", spec$objects),
               check = paste0("doc:usage:", spec$rdname),
               already = paste0("#' ", spec$objects[1]))
  for (o in spec$objects) {
    append_stanza(path, spec$rdname, o, check = paste0("doc:stanza:", o))
  }
}

# ---- Objects withdrawn this cycle --------------------------------------------
# Step 3 drops these by not writing them, which leaves their roxygen behind: the
# family's man page would keep an \alias for an object the package no longer has,
# and `?BLOOD_EPIGEN_ATAC_SEQ_SUM_STATS` would open a page describing data that
# will not load. Both the name in the family's @description listing and the
# object's own stanza have to go.
#
# Driven off the manifest rather than a hand-kept list, and idempotent: a name
# that is already absent is reported, not treated as a missing anchor.

remove_data_doc <- function(path, object) {
  lines <- readLines(path, warn = FALSE)
  quoted <- paste0('"', object, '"')
  n0 <- length(lines)

  # The stanza is the quoted object name preceded by its contiguous roxygen
  # block, plus the blank line that separated it from the stanza above.
  k <- which(lines == quoted)
  if (length(k)) {
    k <- k[1]
    s <- k
    while (s > 1L && grepl("^\\s*#'", lines[s - 1L])) s <- s - 1L
    while (s > 1L && !nzchar(trimws(lines[s - 1L]))) s <- s - 1L
    lines <- lines[-seq.int(s, k)]
  }

  # The name also appears on its own in the family block's tissue listing.
  mention <- which(trimws(lines) == paste0("#' ", object))
  if (length(mention)) lines <- lines[-mention]

  if (length(lines) == n0) return(0L)
  writeLines(lines, path)
  n0 - length(lines)
}

removed <- manifest[manifest$verdict == "REMOVED", , drop = FALSE]
for (i in seq_len(nrow(removed))) {
  obj <- removed$object[i]
  check <- paste0("doc:removed:", obj)
  files <- list.files(file.path(root(removed$destination[i]), "R"),
                      pattern = "[.]R$", full.names = TRUE)
  hits <- vapply(files, function(f) remove_data_doc(f, obj), integer(1))
  if (sum(hits)) {
    record(rpt, "PASS", check,
           sprintf("dropped %d documentation line(s) from %s", sum(hits),
                   paste(basename(files[hits > 0L]), collapse = ", ")))
  } else {
    record(rpt, "PASS", check, "no documentation refers to it")
  }
}

# ---- HUMAN_FEATURE_TO_GENE: regenerate @format from the object itself ---------
# The shipped block described 1,958,568 rows and 11 columns. The object it
# documented had 9 columns; the carried one has 10. Two documented columns are
# gone and one present column was never documented, so this is rewritten from
# what is actually there rather than patched.

COLUMN_DOC <- c(
  assay = "character; the assay (ome). One of \"epigen-atac-seq\", \"epigen-methylcap-seq\", \"metab\", \"prot-clinical\", \"prot-ol\", \"prot-ph\", \"prot-pr\", or \"transcript-rna-seq\".",
  feature_id = "factor; the feature identifier.",
  entrez_gene = "factor; Entrez gene identifier.",
  gene_symbol = "factor; gene symbol.",
  ensembl_gene = "factor; ensembl gene identifier.",
  custom_annotation = "factor; the region of the assigned gene the peak falls in \u2014 one of \"Promoter (<=1kb)\", \"Promoter (1-2kb)\", \"5' UTR\", \"3' UTR\", \"Exon\", \"Intron\", \"Overlaps Gene\", \"Upstream (<5kb)\", \"Downstream (<5kb)\" or \"Distal Intergenic\". \\code{NA} unless assay is \"epigen-atac-seq\" or \"epigen-methylcap-seq\".",
  relationship_to_gene = "numeric; the signed distance in base pairs from the peak to the assigned gene, \\code{0} where the peak overlaps it. \\code{NA} unless assay is \"epigen-atac-seq\" or \"epigen-methylcap-seq\".",
  uniprot = "factor; UniProt identifier.",
  refmet_name = "factor; RefMet metabolite name, from the pinned RefMet snapshot.",
  refmet_id = "factor; RefMet metabolite identifier, from the pinned RefMet snapshot.",
  kegg_id = "factor; KEGG identifier, from the pinned RefMet/KEGG snapshot.",
  flanking_sequence = "factor; flanking sequence. Only applicable if assay is \"prot-ph\"."
)
# The prot-ph / prot-pr per-tissue columns are not columns of this object: confident_site,
# confident_score, ptm_score, redundant_ids, num_peptides, percent_coverage and protein_score
# are all measured per tissue and this table is keyed on (assay, feature_id) with no tissue
# column. They live on \code{*_PROT_PH_QC$feature_metadata} and
# \code{*_PROT_PR_QC$feature_metadata} in MotrpacHumanPreSuspensionData.

h2g_path <- file.path(root("analysis"), "data", "HUMAN_FEATURE_TO_GENE.rda")
if (file.exists(h2g_path)) {
  h2g <- load_rda_object(h2g_path, "HUMAN_FEATURE_TO_GENE")
  undocumented <- setdiff(names(h2g), names(COLUMN_DOC))
  if (length(undocumented)) {
    record(rpt, "FAIL", "doc:HUMAN_FEATURE_TO_GENE",
           paste0("no description available for column(s): ", paste(undocumented, collapse = ", ")))
  } else {
    block <- c(
      sprintf("#' @format A sorted \\code{data.table} with %s rows and %d columns:",
              format(nrow(h2g), big.mark = ","), ncol(h2g)),
      "#'",
      "#' \\describe{",
      unlist(lapply(names(h2g), function(cn) {
        paste0("#'   \\item{", cn, "}{", COLUMN_DOC[[cn]], "}")
      })),
      "#' }",
      "#'",
      "#' @source Built by Stage 1 step 07 of the motrpac-human-presuspension-repro pipeline. The",
      "#'   \\code{refmet_name}, \\code{refmet_id} and \\code{kegg_id} columns come from a",
      "#'   pinned offline RefMet/KEGG snapshot rather than a live Metabolomics Workbench",
      "#'   query, so the mapping does not move with those databases."
    )
    # The rewrite consumes its own anchor, so a re-run has to recognise the
    # already-rewritten form rather than fail looking for the original.
    h2g_r <- file.path(root("analysis"), "R", "HUMAN_FEATURE_TO_GENE.R")
    current <- readLines(h2g_r, warn = FALSE)
    if (any(current == block[1])) {
      record(rpt, "PASS", "doc:HUMAN_FEATURE_TO_GENE",
             "@format already describes the carried object")
    } else {
      replace_block(h2g_r,
                    from = "#' @format A sorted \\code{data.table} with 1,920,618 rows and 10 columns:",
                    to = "#' }", new_lines = block, check = "doc:HUMAN_FEATURE_TO_GENE")
    }
  }
} else {
  record(rpt, "FAIL", "doc:HUMAN_FEATURE_TO_GENE", "object not present in the test package")
}

# ---- No provenance table in either package ------------------------------------
# A copy left in the source checkout is removed so it cannot ship.

for (dest in names(PKG)) {
  prov_path <- file.path(root(dest), "inst", "PROVENANCE.tsv")
  if (file.exists(prov_path)) {
    file.remove(prov_path)
    record(rpt, "PASS", paste0("provenance:", PKG[[dest]]), "removed inst/PROVENANCE.tsv")
  } else {
    record(rpt, "PASS", paste0("provenance:", PKG[[dest]]), "no inst/PROVENANCE.tsv")
  }
}

# ---- METABOLOMICS_CVS: say where it is staged from ---------------------------

mcv <- file.path(root("analysis"), "R", "METABOLOMICS_CVS.R")
if (file.exists(mcv)) {
  lines <- readLines(mcv, warn = FALSE)
  if (!any(grepl("staged verbatim", lines, fixed = TRUE))) {
    i <- grep("^\"METABOLOMICS_CVS\"$", lines)
    if (length(i)) {
      src <- c("#' @source Staged verbatim from the production-bucket metabolite-CV tier;",
               "#'   not regenerated by the motrpac-human-presuspension-repro pipeline.",
               "#'")
      writeLines(append(lines, src, after = i[1] - 1L), mcv)
      record(rpt, "PASS", "doc:METABOLOMICS_CVS", "recorded staged provenance")
    } else {
      record(rpt, "FAIL", "doc:METABOLOMICS_CVS", "object stanza not found")
    }
  } else {
    record(rpt, "PASS", "doc:METABOLOMICS_CVS", "provenance already recorded")
  }
}

# ---- R/run_SCION.R is withdrawn from the analysis package --------------------
# SCION inference is Stage 1 step 16 now, which vendors run_SCION() and its four
# helpers. This step used to patch a guard into R/run_SCION.R, because the four
# prot-ph / prot-pr QC objects drop `qc_imputed` and the function assigned
#   x[[tissue]][[ome]][["qc_norm"]] <- x[[tissue]][[ome]][["qc_imputed"]]
# — assigning NULL to a list element removes it, so with the component gone that
# line deleted qc_norm rather than replacing it. Step 16 reads the imputed matrix
# from step 06's freeze output instead, so there is nothing left to guard.
#
# The file's absence is asserted rather than assumed: a source checkout that still
# carries it would ship four exports whose helpers reach into the Data package
# through string literals, an undeclared runtime cycle.

scion <- file.path(root("analysis"), "R", "run_SCION.R")

if (!file.exists(scion)) {
  record(rpt, "PASS", "code:run_SCION", "withdrawn from the package; inference is Stage 1 step 16")
} else {
  record(rpt, "FAIL", "code:run_SCION",
         paste("R/run_SCION.R is back in the analysis package — it was withdrawn in favour of",
               "Stage 1 step 16; delete it, its two man/ pages and its four NAMESPACE exports"))
}

# ---- plot_single_feature(), onto the split clinical-chemistry objects --------
# v1.3 carried clinical chemistry as one assay and step 2 withdraws that pair, so
# the three shipped call sites that read them have to move to the split. Left
# alone this is not a stale reference but a hard failure: the objects are gone,
# and `Pkg::MISSING` errors at call time.
#
# The rewrite preserves what the function renders today. The split changes the
# assay labels as well as the object names — the combined objects carried
# `assay = "clinical-chemistry"`, where the halves carry `metab` (with
# `platform = "metab-t-clinical"`) and `prot-clinical` — so two things follow:
#
#  * the DA half needs the same metab-platform normalisation the function already
#    applies to its main DA frame, or `assay` disagrees between the DA and the
#    summary statistics and the full_join on (tissue, assay, Timepoint,
#    randomGroupCode, feature_id) matches nothing.
#  * the "Clin. Chem." facet label was keyed on the old assay string, so it has
#    to recognise the two new ones or clinical panels come out labelled NA.
#
# The split is an exact partition — 198 + 99 = 297 DA rows, 114 + 57 = 171
# summary-stat rows, the same totals the combined objects held — so binding the
# halves reproduces the old input row for row.

psf <- file.path(root("analysis"), "R", "plot_single_feature.R")
psf_src <- if (is.null(opt$`analysis-pkg`)) NULL else
  file.path(opt$`analysis-pkg`, "R", "plot_single_feature.R")
PSF_MARK <- "the v2.0 split of clinical chemistry"

if (!file.exists(psf)) {
  record(rpt, "WARN", "code:plot_single_feature", "plot_single_feature.R not found; nothing to adapt")
} else {
  if (!is.null(psf_src) && file.exists(psf_src)) file.copy(psf_src, psf, overwrite = TRUE)
  lines <- readLines(psf, warn = FALSE)

  if (any(grepl(PSF_MARK, lines, fixed = TRUE))) {
    record(rpt, "PASS", "code:plot_single_feature", "already on the split objects — no change")
  } else {
    # Each edit names the exact line it replaces and fails if it is not there,
    # so a silent no-match cannot pass for a completed rewrite.
    psf_edits <- list(
      list(what = "da",
           pat = "^(\\s*)clin_chem_da_rows = MotrpacHumanPreSuspensionAnalysis::CLIN_CHEMISTRY_DA %>%$",
           new = function(ind) c(
             paste0(ind, "# ", PSF_MARK, " replaces CLIN_CHEMISTRY_DA with one object per"),
             paste0(ind, "# assay. The mutate is the same metab-platform normalisation applied to the"),
             paste0(ind, "# main DA frame above, and it is what makes `assay` agree with the summary"),
             paste0(ind, "# statistics these rows are later joined to."),
             paste0(ind, "clin_chem_da_rows = dplyr::bind_rows("),
             paste0(ind, "    MotrpacHumanPreSuspensionAnalysis::BLOOD_METAB_T_CLINICAL_DA,"),
             paste0(ind, "    MotrpacHumanPreSuspensionAnalysis::BLOOD_PROT_CLINICAL_DA"),
             paste0(ind, "  ) %>%"),
             paste0(ind, "  dplyr::mutate(assay = ifelse(assay == \"metab\", as.character(platform), as.character(assay))) %>%"))),
      list(what = "sum",
           pat = "^(\\s*)clin_chem_sum = MotrpacHumanPreSuspensionAnalysis::BLOOD_CLINICAL_CHEMISTRY_SUM_STATS %>%$",
           new = function(ind) c(
             paste0(ind, "clin_chem_sum = dplyr::bind_rows("),
             paste0(ind, "    MotrpacHumanPreSuspensionAnalysis::BLOOD_METAB_T_CLINICAL_SUM_STATS,"),
             paste0(ind, "    MotrpacHumanPreSuspensionAnalysis::BLOOD_PROT_CLINICAL_SUM_STATS"),
             paste0(ind, "  ) %>%"))),
      list(what = "ids",
           pat = "^(\\s*)ids = MotrpacHumanPreSuspensionAnalysis::BLOOD_CLINICAL_CHEMISTRY_SUM_STATS\\$feature_id$",
           new = function(ind) c(
             paste0(ind, "ids = c(MotrpacHumanPreSuspensionAnalysis::BLOOD_METAB_T_CLINICAL_SUM_STATS$feature_id,"),
             paste0(ind, "        MotrpacHumanPreSuspensionAnalysis::BLOOD_PROT_CLINICAL_SUM_STATS$feature_id)"))),
      list(what = "label",
           pat = "^(\\s*)assay_short_text = dplyr::if_else\\(assay == \"clinical-chemistry\", \"Clin. Chem.\", assay_short_text\\),$",
           new = function(ind) c(
             paste0(ind, "assay_short_text = dplyr::if_else("),
             paste0(ind, "  assay %in% c(\"clinical-chemistry\", \"metab-t-clinical\", \"prot-clinical\"),"),
             paste0(ind, "  \"Clin. Chem.\", assay_short_text),"))),
      list(what = "doc",
           pat = "^(#')  *\\\\code\\{feature_id\\} column of \\\\code\\{BLOOD_CLINICAL_CHEMISTRY_SUM_STATS\\}\\.$",
           new = function(ind) paste0(ind, " \\code{feature_id} columns of \\code{BLOOD_METAB_T_CLINICAL_SUM_STATS} and \\code{BLOOD_PROT_CLINICAL_SUM_STATS}."))
    )

    for (e in psf_edits) {
      k <- grep(e$pat, lines)
      if (!length(k)) {
        record(rpt, "FAIL", paste0("code:plot_single_feature:", e$what),
               paste0("anchor not found: ", e$pat))
        next
      }
      k <- k[1]
      ind <- sub(e$pat, "\\1", lines[k])
      lines <- append(lines[-k], e$new(ind), after = k - 1L)
      record(rpt, "PASS", paste0("code:plot_single_feature:", e$what),
             "moved onto the split clinical-chemistry objects")
    }
    writeLines(lines, psf)
  }
}

# The load_qc() roxygen still describes qc_imputed as a element users can expect
# for prot-pr/ph. No object ships it now.
lq <- file.path(root("data"), "R", "load_qc.R")
if (file.exists(lq)) {
  lines <- readLines(lq, warn = FALSE)
  hits <- grep("gone through multiple imputation using \\\\code\\{run_mice\\}. Only for prot-pr/ph", lines)
  if (!length(hits)) {
    record(rpt, "PASS", "doc:load_qc:qc_imputed", "no stale qc_imputed description found")
  } else {
    lines[hits] <- sub("Only for prot-pr/ph",
                       "No longer shipped as of this release; present only in archived objects",
                       lines[hits])
    writeLines(lines, lq)
    record(rpt, "PASS", "doc:load_qc:qc_imputed",
           sprintf("updated %d description(s) of the dropped qc_imputed element", length(hits)))
  }
}

# ---- Every R file this step touched must still parse -------------------------
# Editing R source with line surgery can produce a file that is syntactically
# valid but structurally wrong — inserting under an unbraced `if` rebinds its
# body. Parsing catches the first kind; the block rewrite above is what handles
# the second. Both are worth having, and parsing has to happen before roxygen,
# which would otherwise report the same problem far less clearly.

for (dest in names(PKG)) {
  bad <- character()
  for (f in list.files(file.path(root(dest), "R"), pattern = "[.]R$", full.names = TRUE)) {
    ok <- tryCatch({ parse(f); TRUE }, error = function(e) FALSE)
    if (!ok) bad <- c(bad, basename(f))
  }
  if (length(bad)) {
    record(rpt, "FAIL", paste0("parse:", PKG[[dest]]),
           paste0("R file(s) do not parse: ", paste(bad, collapse = ", ")))
  } else {
    record(rpt, "PASS", paste0("parse:", PKG[[dest]]), "every R file parses")
  }
}

# ---- Re-roxygenise ------------------------------------------------------------
# man/load_differential_analysis.Rd still documents CI.L and CI.R; the roxygen
# source that generated it no longer does. Regenerating is what drops them.
#
# roxygenise() MUST run with the working directory set to the package root. The
# 76 cln_* stubs and pheno document their variables through
# `@includeRmd rmd/<object>.Rmd`, and every one of those Rmd opens with
# ```{r child = 'aaa.Rmd'}```. knitr resolves a child path against the WORKING
# DIRECTORY, not against the package or the parent document, so rendering from
# anywhere else fails to find aaa.Rmd, the render yields nothing, and each page
# collapses from ~780 lines to the 21 that roxygen can produce without it — the
# whole per-variable data dictionary, silently. roxygen2 reports no error and the
# man-page COUNT is unchanged, so nothing downstream notices. Only the shrink
# guard below would.

#' Line count of every man page, keyed by file name
man_page_lines <- function(p) {
  files <- list.files(file.path(p, "man"), pattern = "\\.Rd$", full.names = TRUE)
  out <- vapply(files, function(f) length(readLines(f, warn = FALSE)), integer(1))
  names(out) <- basename(files)
  return(out)
}

if (!requireNamespace("roxygen2", quietly = TRUE)) {
  record(rpt, "FAIL", "roxygen", "roxygen2 is not installed")
} else {
  for (dest in names(PKG)) {
    p <- root(dest)
    before <- man_page_lines(p)

    res <- tryCatch({
      # setwd rather than withr, which neither package depends on.
      old_wd <- setwd(p)
      on.exit(setwd(old_wd), add = TRUE)
      suppressMessages(roxygen2::roxygenise(".", roclets = c("namespace", "rd")))
      setwd(old_wd)
      TRUE
    }, error = function(e) { record(rpt, "FAIL", paste0("roxygen:", PKG[[dest]]), conditionMessage(e)); FALSE })

    if (isTRUE(res)) {
      n <- length(list.files(file.path(p, "man"), pattern = "\\.Rd$"))
      record(rpt, "PASS", paste0("roxygen:", PKG[[dest]]), sprintf("%d man page(s)", n))

      # A regenerated page is allowed to lose a little — a dropped @format entry
      # is a line or two. Losing most of a page means its generated content did
      # not render, which is a silent failure by construction: the file is still
      # valid Rd and R CMD check reads it happily.
      after  <- man_page_lines(p)
      shared <- base::intersect(names(before), names(after))
      shrunk <- shared[after[shared] < 0.5 * before[shared] & before[shared] >= 40]
      if (length(shrunk)) {
        record(rpt, "FAIL", paste0("roxygen:shrink:", PKG[[dest]]),
               sprintf("%d man page(s) lost more than half their content, so generated documentation did not render: %s",
                       length(shrunk),
                       paste(sprintf("%s (%d->%d)", head(shrunk, 5), before[head(shrunk, 5)],
                                     after[head(shrunk, 5)]), collapse = ", ")))
      } else {
        record(rpt, "PASS", paste0("roxygen:shrink:", PKG[[dest]]),
               sprintf("no man page lost more than half its content (%d compared)", length(shared)))
      }
    }
  }
}

# ---- Verify: every data object has a man page, and vice versa ----------------

for (dest in names(PKG)) {
  p <- root(dest)
  objs <- sub("\\.rda$", "", list.files(file.path(p, "data"), pattern = "\\.rda$"))
  rd <- list.files(file.path(p, "man"), pattern = "\\.Rd$", full.names = TRUE)
  aliases <- unlist(lapply(rd, function(f) {
    l <- readLines(f, warn = FALSE)
    sub("^\\\\alias\\{(.*)\\}$", "\\1", l[grepl("^\\\\alias\\{", l)])
  }))
  undoc <- setdiff(objs, aliases)
  if (length(undoc)) {
    record(rpt, "FAIL", paste0("doc:coverage:", PKG[[dest]]),
           sprintf("%d undocumented data object(s): %s", length(undoc),
                   paste(undoc, collapse = ", ")))
  } else {
    record(rpt, "PASS", paste0("doc:coverage:", PKG[[dest]]),
           sprintf("all %d data object(s) documented", length(objs)))
  }

  # The other direction, but only for the objects this cycle withdrew. Aliases
  # with no object behind them are normal — every exported function has one — so
  # a blanket check here would be noise.
  gone <- manifest$object[manifest$destination == dest & manifest$verdict == "REMOVED"]
  dangling <- intersect(gone, aliases)
  if (length(dangling)) {
    record(rpt, "FAIL", paste0("doc:withdrawn:", PKG[[dest]]),
           sprintf("man page(s) still alias withdrawn object(s): %s",
                   paste(dangling, collapse = ", ")))
  } else if (length(gone)) {
    record(rpt, "PASS", paste0("doc:withdrawn:", PKG[[dest]]),
           sprintf("no man page aliases the %d withdrawn object(s)", length(gone)))
  }
}

# CI.L and CI.R must be gone from the regenerated man pages: they are no longer
# columns of any DA object. CI.L_calculated / CI.R_calculated ARE columns and may be
# documented, so the pattern excludes them — a plain "CI.L" substring match would read
# the corrected pair as the removed one and fail every page that documents it.
rd <- list.files(file.path(root("analysis"), "man"), pattern = "\\.Rd$", full.names = TRUE)
stale <- rd[vapply(rd, function(f) any(grepl("CI\\.[LR](?!_calculated)",
                                             readLines(f, warn = FALSE), perl = TRUE)),
                   logical(1))]
if (length(stale)) {
  record(rpt, "FAIL", "doc:CI-columns",
         paste0("man page(s) still document the removed CI.L/CI.R columns: ",
                paste(basename(stale), collapse = ", ")))
} else {
  record(rpt, "PASS", "doc:CI-columns", "no man page documents the removed CI.L/CI.R columns")
}

write_tsv(report_frame(rpt), file.path(opt$out, "document_report.tsv"))
message(sprintf("\n%d check(s): %d FAIL, %d WARN", length(rpt$rows), rpt$fails, rpt$warns))
if (rpt$fails > 0L) quit(status = 1L)
