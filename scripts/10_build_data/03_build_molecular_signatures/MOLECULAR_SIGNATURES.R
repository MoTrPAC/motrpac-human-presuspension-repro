#!/usr/bin/env Rscript
# Stage 1 data object: MOLECULAR_SIGNATURES  (deps: gmt_files)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/MOLECULAR_SIGNATURES.R.
#
# Local-.rda build: instead of the package's data-raw/gmt_files/, it reads the 7
# GMTs already built by preflight (scripts/00_preflight/output/), and instead of
# usethis::use_data() it save()s a LOCAL .rda into scripts/10_build_data/data/.
# There is no in-house library()/:: dependency to replace here — the only object
# input is the gmt_files node, consumed as those local .gmt.gz.
#
# NOTE: this data/MOLECULAR_SIGNATURES.rda is staging — it will eventually be
# organized as the MOLECULAR_SIGNATURES data object in the Analysis package.
suppressWarnings(suppressMessages(library(TMSig)))  # readGMT

# Root from the PRECOVID_ROOT config env (config/pipeline.env); fall back to an
# upward search for config/pipeline.env so standalone runs work.
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
}
gmt_dir <- file.path(ROOT, "scripts", "00_preflight", "output")
out_dir <- file.path(ROOT, "scripts", "10_build_data", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gmt_files <- c(
  "C2"         = "c2.all.v2023.2.Hs.symbols.gmt.gz",                       # MSigDB C2
  "GO"         = "c5.go.v2023.2.Hs.symbols.gmt.gz",                        # MSigDB C5/GO
  "MITOCARTA"  = "mitocarta3.0.symbols.gmt.gz",                            # MitoCarta 3.0
  "PSP"        = "phosphositeplus.v6.7.1.1.flanking.gmt.gz",               # kinase sets
  "REFMET"     = "metabolomics.workbench.refmet.2024.08.07.metabolites.gmt.gz",
  "PTMSIGDB"   = "ptmsigdb.v2.0.flanking.gmt.gz",                          # PTMSigDB
  "CELLMARKER" = "cellmarker.v2024.symbols.gmt.gz"                         # CellMarker
)
gmt_files <- structure(file.path(gmt_dir, gmt_files), names = names(gmt_files))
missing <- gmt_files[!file.exists(gmt_files)]
if (length(missing))
  stop("missing GMT input(s) — build preflight data objects first:\n  ",
       paste(missing, collapse = "\n  "))

# Nested list of molecular signatures
MOLECULAR_SIGNATURES <- lapply(gmt_files, readGMT)

# Separate Gene Ontology databases
GO <- MOLECULAR_SIGNATURES[["GO"]]
ont <- sub("^(GO[^_]+)_.*", "\\1", names(GO))
GO <- split(do.call(list, GO), f = ont)

C2 <- MOLECULAR_SIGNATURES[["C2"]]
databases <- c("BIOCARTA", "KEGG_MEDICUS", "PID", "REACTOME", "WP")
keep <- grepl(paste(paste0("^", databases, "_"), collapse = "|"), names(C2))
C2 <- C2[keep]
group <- sub("^([^_]+).*", "\\1", names(C2))
group[group == "KEGG"] <- "KEGG_MEDICUS"
C2 <- split(do.call(list, C2), f = group)

MOLECULAR_SIGNATURES[c("C2", "GO")] <- NULL
MOLECULAR_SIGNATURES <- c(C2, GO, MOLECULAR_SIGNATURES)

save(MOLECULAR_SIGNATURES, file = file.path(out_dir, "MOLECULAR_SIGNATURES.rda"),
     compress = TRUE, version = 3)
message("MOLECULAR_SIGNATURES: ", length(MOLECULAR_SIGNATURES), " databases -> ",
        file.path(out_dir, "MOLECULAR_SIGNATURES.rda"))
