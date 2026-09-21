# Run environment

**Generated automatically — do not edit by hand.** Regenerate on demand with `make env`; it is not
part of the pipeline DAG and does not re-run on every build. Refresh it when the environment changes
or before committing regenerated outputs: it is the record of *which* machine, *when*, and *with
which package versions* those outputs were produced.

Companion machine-readable files live in `docs/environment/`:
`package_versions.tsv` (full package inventory) and `session_info.txt` (raw R `sessionInfo()`).

## Run metadata

| Field | Value |
|---|---|
| Date of run | **2026-07-30 09:46:30 CDT (-0500)** |
| Unix timestamp | 1785422790 |
| User | christopherjin |
| Host | Christophers-MacBook-Pro-5.local |
| Pipeline root | `/Users/christopherjin/Documents/GitHub/motrpac-human-presuspension-repro` |
| Data versions | current `v1.3` → new `v1.4` |
| Locale | en_US.UTF-8 |

## Operating system

| Field | Value |
|---|---|
| OS | macOS 26.5.2 (build 25F84) |
| Kernel | Darwin 25.5.0 |
| Architecture | arm64 |
| CPU | Apple M2 Max |
| Cores | 12 |
| Memory | 32 GB |

## Toolchain

| Tool | Version | Path |
|---|---|---|
| R | R version 4.4.1 (2024-06-14) | `/usr/local/bin/Rscript` |
| bash | 3.2.57(1)-release | `/bin/bash` |
| make | GNU Make 3.81 | `/usr/bin/make` |
| git | git version 2.39.1 | `/opt/homebrew/bin/git` |
| gsutil | gsutil version: 5.24 | `/Users/christopherjin/google-cloud-sdk/bin/gsutil` |
| gcloud | Google Cloud SDK 435.0.1 | `/Users/christopherjin/google-cloud-sdk/bin/gcloud` |
| curl | curl 8.7.1 (x86_64-apple-darwin25.0) libcurl/8.7.1 (SecureTransport) LibreSSL/3.3.6 zlib/1.2.12 nghttp2/1.68.1 | `/usr/bin/curl` |

## Repository state

Every repo the pipeline reads or writes, at the commit used for this run.

| Repo | Package version | Git state |
|---|---|---|
| `motrpac-human-presuspension-repro` | — | `main` @ `708cd0e` (**13 uncommitted file(s)**) |
| `MotrpacHumanPreSuspensionData` | 0.0.1.102 | `updating_feature_metadata` @ `4aa922f` (clean) |
| `MotrpacHumanPreSuspensionAnalysis` | 0.2.4 | `clarifying_dependencies` @ `2dbb25f` (clean) |
| `MotrpacPreSuspensionAcute` | — | `chris_revisions` @ `6b41684` (**16 uncommitted file(s)**) |

## R package versions

436 packages resolved: 118 declared directly by the three in-scope repos (or checked by
preflight), the rest pulled in as recursive `Depends`/`Imports`/`LinkingTo`. Versions are what is
installed on this machine right now. Full table, including install library paths:
[`docs/environment/package_versions.tsv`](environment/package_versions.tsv).

### Directly declared (118)

| Package | Version | Source | Built under |
|---|---|---|---|
| `AnnotationDbi` | 1.68.0 | Bioconductor 3.20 | 4.4.1 |
| `annotatr` | 1.32.0 | Bioconductor 3.20 | 4.4.1 |
| `base` | 4.4.1 | base | 4.4.1 |
| `Biobase` | 2.72.0 | Bioconductor | 4.4.1 |
| `BiocManager` | 1.30.26 | CRAN | 4.4.1 |
| `BiocParallel` | 1.38.0 | Bioconductor | 4.4.0 |
| `biomaRt` | 2.60.1 | Bioconductor | 4.4.0 |
| `bubbleHeatmap` | 0.1.1 | CRAN | 4.4.1 |
| `Cairo` | 1.6.2 | CRAN | 4.4.1 |
| `car` | 3.1.5 | CRAN | 4.4.3 |
| `ChIPseeker` | 1.42.1 | Bioconductor 3.20 | 4.4.2 |
| `circlize` | 0.4.17 | CRAN | 4.4.3 |
| `clusterProfiler` | 4.12.6 | Bioconductor | 4.4.1 |
| `cmapR` | 1.16.0 | Bioconductor | 4.4.0 |
| `ComplexHeatmap` | 2.20.0 | Bioconductor | 4.4.0 |
| `ComplexUpset` | 1.3.6 | github:krassowski/complex-upset@79fa000 | 4.4.1 |
| `conflicted` | 1.2.0 | CRAN | 4.4.0 |
| `cowplot` | 1.2.0 | CRAN | 4.4.1 |
| `curl` | 7.0.0 | CRAN | 4.4.1 |
| `data.table` | 1.18.2.1 | CRAN | 4.4.3 |
| `devtools` | 2.4.5 | CRAN | 4.4.0 |
| `doParallel` | 1.0.17 | CRAN | 4.4.0 |
| `dplyr` | 1.2.0 | CRAN | 4.4.3 |
| `edgeR` | 4.2.2 | Bioconductor | 4.4.1 |
| `ensembldb` | 2.30.0 | Bioconductor 3.20 | 4.4.1 |
| `fgsea` | 1.30.0 | Bioconductor | 4.4.0 |
| `forcats` | 1.0.1 | CRAN | 4.4.1 |
| `foreach` | 1.5.2 | CRAN | 4.4.0 |
| `futile.logger` | 1.4.9 | CRAN | 4.4.3 |
| `gdsfmt` | 1.42.1 | Bioconductor 3.20 | 4.4.2 |
| `GenomeInfoDb` | 1.40.1 | Bioconductor | 4.4.0 |
| `GenomicRanges` | 1.58.0 | Bioconductor 3.20 | 4.4.1 |
| `geomtextpath` | 0.2.0 | CRAN | 4.4.1 |
| `ggbeeswarm` | 0.7.2 | CRAN | 4.4.0 |
| `ggh4x` | 0.3.1 | CRAN | 4.4.1 |
| `ggnewscale` | 0.5.0 | CRAN | 4.4.0 |
| `ggpattern` | 1.3.1 | CRAN | 4.4.3 |
| `ggplot2` | 4.0.2 | CRAN | 4.4.3 |
| `ggpubr` | 0.6.3 | CRAN | 4.4.3 |
| `ggrastr` | 1.0.2 | CRAN | 4.4.0 |
| `ggrepel` | 0.9.8 | CRAN | 4.4.3 |
| `ggsankey` | 0.0.99999 | github:davidsjoberg/ggsankey@b675d0d | 4.4.1 |
| `ggtext` | 0.1.2 | CRAN | 4.4.0 |
| `ggVennDiagram` | 1.5.2 | CRAN | 4.4.0 |
| `glue` | 1.8.1 | CRAN | 4.4.1 |
| `gplots` | 3.3.0 | CRAN | 4.4.3 |
| `gprofiler2` | 0.2.3 | CRAN | 4.4.0 |
| `grDevices` | 4.4.1 | base | 4.4.1 |
| `grid` | 4.4.1 | base | 4.4.1 |
| `gridExtra` | 2.3 | CRAN | 4.4.1 |
| `gridtext` | 0.1.5 | CRAN | 4.4.0 |
| `GSEABase` | 1.66.0 | Bioconductor | 4.4.0 |
| `here` | 1.0.2 | CRAN | 4.4.1 |
| `igraph` | 2.2.1 | CRAN | 4.4.1 |
| `IRanges` | 2.40.1 | Bioconductor 3.20 | 4.4.2 |
| `jsonlite` | 2.0.0 | CRAN | 4.4.1 |
| `KEGGREST` | 1.44.1 | Bioconductor | 4.4.0 |
| `knitr` | 1.51 | CRAN | 4.4.3 |
| `latex2exp` | 0.9.8 | CRAN | 4.4.3 |
| `limma` | 3.62.2 | Bioconductor | 4.4.2 |
| `lme4` | 1.1.38 | CRAN | 4.4.1 |
| `lmerTest` | 3.2.1 | CRAN | 4.4.3 |
| `magrittr` | 2.0.5 | CRAN | 4.4.1 |
| `matrixStats` | 1.5.0 | CRAN | 4.4.1 |
| `Mfuzz` | 2.72.0 | Bioconductor | 4.4.1 |
| `mice` | 3.17.0 | CRAN | 4.4.1 |
| `MotrpacBicQC` | 1.7.0 | github:MoTrPAC/MotrpacBicQC@d873a93 | 4.4.1 |
| `MotrpacHumanPreSuspensionAnalysis` | 0.2.4 | local source install | 4.4.1 |
| `MotrpacHumanPreSuspensionData` | 0.0.1.102 | local source install | 4.4.1 |
| `MotrpacRatTraining6mo` | 1.6.5 | github:MoTrPAC/MotrpacRatTraining6mo@0f38357 | 4.4.1 |
| `msigdbr` | 25.1.1 | CRAN | 4.4.1 |
| `openxlsx` | 4.2.8 | CRAN | 4.4.1 |
| `org.Hs.eg.db` | 3.19.1 | Bioconductor | 4.4.1 |
| `pak` | 0.8.0 | CRAN | 4.4.1 |
| `parallel` | 4.4.1 | base | 4.4.1 |
| `patchwork` | 1.3.2 | CRAN | 4.4.1 |
| `pheatmap` | 1.0.12 | CRAN | 4.4.0 |
| `pkgdown` | 2.1.1 | CRAN | 4.4.1 |
| `PLIER` | 0.99.0 | github:wgmao/PLIER@fe4e9b2 | 4.4.1 |
| `psych` | 2.4.6.26 | CRAN | 4.4.0 |
| `purrr` | 1.2.1 | CRAN | 4.4.3 |
| `R.utils` | 2.12.3 | CRAN | 4.4.0 |
| `randomForest` | 4.7.1.2 | CRAN | 4.4.1 |
| `RColorBrewer` | 1.1.3 | CRAN | 4.4.0 |
| `readr` | 2.2.0 | CRAN | 4.4.3 |
| `readxl` | 1.4.3 | CRAN | 4.4.0 |
| `reformulas` | 0.4.4 | CRAN | 4.4.3 |
| `remotes` | 2.5.0 | CRAN | 4.4.0 |
| `renv` | 1.0.11 | CRAN | 4.4.1 |
| `reshape2` | 1.4.5 | CRAN | 4.4.1 |
| `rjson` | 0.2.23 | CRAN | 4.4.1 |
| `rlang` | 1.3.0 | CRAN | 4.4.1 |
| `rmarkdown` | 2.30 | CRAN | 4.4.1 |
| `Rmisc` | 1.5.1 | CRAN | 4.4.0 |
| `roxygen2` | 8.0.0 | CRAN | 4.4.1 |
| `Rtsne` | 0.17 | CRAN | 4.4.0 |
| `scales` | 1.4.0 | CRAN | 4.4.1 |
| `SNPRelate` | 1.40.0 | Bioconductor 3.20 | 4.4.1 |
| `stats` | 4.4.1 | base | 4.4.1 |
| `stringi` | 1.8.7 | CRAN | 4.4.1 |
| `stringr` | 1.6.0 | CRAN | 4.4.1 |
| `table.glue` | 0.0.5 | CRAN | 4.4.1 |
| `testthat` | 3.3.2 | CRAN | 4.4.3 |
| `tibble` | 3.3.1 | CRAN | 4.4.3 |
| `tidyr` | 1.3.2 | CRAN | 4.4.3 |
| `tidyverse` | 2.0.0 | CRAN | 4.4.0 |
| `TMSig` | 1.6.0 | Bioconductor | 4.4.1 |
| `tools` | 4.4.1 | base | 4.4.1 |
| `TxDb.Hsapiens.UCSC.hg19.knownGene` | 3.2.2 | Bioconductor 3.20 | 4.4.1 |
| `txdbmaker` | 1.2.1 | Bioconductor 3.20 | 4.4.2 |
| `umap` | 0.2.10.0 | CRAN | 4.4.0 |
| `UpSetR` | 1.4.0 | CRAN | 4.4.0 |
| `usethis` | 3.1.0 | CRAN | 4.4.1 |
| `utils` | 4.4.1 | base | 4.4.1 |
| `variancePartition` | 1.38.1 | Bioconductor | 4.4.1 |
| `VennDiagram` | 1.7.3 | CRAN | 4.4.0 |
| `WGCNA` | 1.73 | CRAN | 4.4.1 |
| `writexl` | 1.5.4 | CRAN | 4.4.1 |

<details>
<summary><b>All 436 packages, including transitive dependencies</b></summary>

| Package | Version | Source | Built under |
|---|---|---|---|
| `abind` | 1.4.8 | CRAN | 4.4.1 |
| `admisc` | 0.37 | CRAN | 4.4.1 |
| `annotate` | 1.82.0 | Bioconductor | 4.4.0 |
| `AnnotationDbi` | 1.68.0 | Bioconductor 3.20 | 4.4.1 |
| `AnnotationFilter` | 1.30.0 | Bioconductor 3.20 | 4.4.1 |
| `AnnotationHub` | 3.14.0 | Bioconductor 3.20 | 4.4.1 |
| `annotatr` | 1.32.0 | Bioconductor 3.20 | 4.4.1 |
| `aod` | 1.3.3 | CRAN | 4.4.0 |
| `ape` | 5.8 | CRAN | 4.4.0 |
| `aplot` | 0.2.3 | CRAN | 4.4.0 |
| `askpass` | 1.2.1 | CRAN | 4.4.1 |
| `assertthat` | 0.2.1 | CRAN | 4.4.0 |
| `babelgene` | 22.9 | CRAN | 4.4.0 |
| `backports` | 1.5.0 | CRAN | 4.4.1 |
| `base` | 4.4.1 | base | 4.4.1 |
| `base64enc` | 0.1.3 | CRAN | 4.4.0 |
| `beeswarm` | 0.4.0 | CRAN | 4.4.1 |
| `BH` | 1.90.0.1 | CRAN | 4.4.3 |
| `Biobase` | 2.72.0 | Bioconductor | 4.4.1 |
| `BiocFileCache` | 2.12.0 | Bioconductor | 4.4.0 |
| `BiocGenerics` | 0.52.0 | Bioconductor 3.20 | 4.4.1 |
| `BiocIO` | 1.14.0 | Bioconductor | 4.4.0 |
| `BiocManager` | 1.30.26 | CRAN | 4.4.1 |
| `BiocParallel` | 1.38.0 | Bioconductor | 4.4.0 |
| `BiocVersion` | 3.19.1 | Bioconductor | 4.4.0 |
| `biomaRt` | 2.60.1 | Bioconductor | 4.4.0 |
| `Biostrings` | 2.72.1 | Bioconductor | 4.4.0 |
| `bit` | 4.6.0 | CRAN | 4.4.1 |
| `bit64` | 4.6.0.1 | CRAN | 4.4.1 |
| `bitops` | 1.0.9 | CRAN | 4.4.1 |
| `blob` | 1.3.0 | CRAN | 4.4.3 |
| `boot` | 1.3.30 | CRAN | 4.4.1 |
| `brew` | 1.0.10 | CRAN | 4.4.0 |
| `brio` | 1.1.5 | CRAN | 4.4.0 |
| `broom` | 1.0.12 | CRAN | 4.4.3 |
| `BSgenome` | 1.74.0 | Bioconductor 3.20 | 4.4.1 |
| `bslib` | 0.10.0 | CRAN | 4.4.3 |
| `bubbleHeatmap` | 0.1.1 | CRAN | 4.4.1 |
| `cachem` | 1.1.0 | CRAN | 4.4.0 |
| `Cairo` | 1.6.2 | CRAN | 4.4.1 |
| `callr` | 3.7.6 | CRAN | 4.4.0 |
| `car` | 3.1.5 | CRAN | 4.4.3 |
| `carData` | 3.0.6 | CRAN | 4.4.3 |
| `caTools` | 1.18.3 | CRAN | 4.4.1 |
| `cellranger` | 1.1.0 | CRAN | 4.4.0 |
| `checkmate` | 2.3.2 | CRAN | 4.4.0 |
| `ChIPseeker` | 1.42.1 | Bioconductor 3.20 | 4.4.2 |
| `circlize` | 0.4.17 | CRAN | 4.4.3 |
| `class` | 7.3.22 | CRAN | 4.4.1 |
| `classInt` | 0.4.11 | CRAN | 4.4.1 |
| `cli` | 3.6.6 | CRAN | 4.4.1 |
| `clipr` | 0.8.0 | CRAN | 4.4.1 |
| `clue` | 0.3.67 | CRAN | 4.4.3 |
| `cluster` | 2.1.6 | CRAN | 4.4.1 |
| `clusterProfiler` | 4.12.6 | Bioconductor | 4.4.1 |
| `cmapR` | 1.16.0 | Bioconductor | 4.4.0 |
| `codetools` | 0.2.20 | CRAN | 4.4.1 |
| `colorspace` | 2.1.2 | CRAN | 4.4.1 |
| `commonmark` | 2.0.0 | CRAN | 4.4.1 |
| `compiler` | 4.4.1 | base | 4.4.1 |
| `ComplexHeatmap` | 2.20.0 | Bioconductor | 4.4.0 |
| `ComplexUpset` | 1.3.6 | github:krassowski/complex-upset@79fa000 | 4.4.1 |
| `conflicted` | 1.2.0 | CRAN | 4.4.0 |
| `corpcor` | 1.6.10 | CRAN | 4.4.0 |
| `corrplot` | 0.95 | CRAN | 4.4.1 |
| `cowplot` | 1.2.0 | CRAN | 4.4.1 |
| `cpp11` | 0.5.5 | CRAN | 4.4.1 |
| `crayon` | 1.5.3 | CRAN | 4.4.0 |
| `credentials` | 2.0.2 | CRAN | 4.4.1 |
| `crosstalk` | 1.2.1 | CRAN | 4.4.0 |
| `curl` | 7.0.0 | CRAN | 4.4.1 |
| `cytolib` | 2.16.0 | Bioconductor | 4.4.0 |
| `data.table` | 1.18.2.1 | CRAN | 4.4.3 |
| `DBI` | 1.3.0 | CRAN | 4.4.3 |
| `dbplyr` | 2.5.0 | CRAN | 4.4.0 |
| `DelayedArray` | 0.32.0 | Bioconductor 3.20 | 4.4.1 |
| `DEoptimR` | 1.1.3.1 | CRAN | 4.4.1 |
| `Deriv` | 4.2.0 | CRAN | 4.4.1 |
| `desc` | 1.4.3 | CRAN | 4.4.0 |
| `devtools` | 2.4.5 | CRAN | 4.4.0 |
| `diffobj` | 0.3.6 | CRAN | 4.4.1 |
| `digest` | 0.6.39 | CRAN | 4.4.3 |
| `doBy` | 4.7.1 | CRAN | 4.4.3 |
| `doParallel` | 1.0.17 | CRAN | 4.4.0 |
| `DOSE` | 3.30.5 | Bioconductor | 4.4.1 |
| `downlit` | 0.4.4 | CRAN | 4.4.0 |
| `downloader` | 0.4 | CRAN | 4.4.0 |
| `dplyr` | 1.2.0 | CRAN | 4.4.3 |
| `dtplyr` | 1.3.1 | CRAN | 4.4.0 |
| `dynamicTreeCut` | 1.63.1 | CRAN | 4.4.1 |
| `DynDoc` | 1.82.0 | Bioconductor | 4.4.0 |
| `e1071` | 1.7.16 | CRAN | 4.4.1 |
| `edgeR` | 4.2.2 | Bioconductor | 4.4.1 |
| `ellipsis` | 0.3.2 | CRAN | 4.4.0 |
| `enrichplot` | 1.24.4 | Bioconductor | 4.4.1 |
| `ensembldb` | 2.30.0 | Bioconductor 3.20 | 4.4.1 |
| `EnvStats` | 3.1.0 | CRAN | 4.4.1 |
| `evaluate` | 1.0.5 | CRAN | 4.4.1 |
| `fANCOVA` | 0.6.1 | CRAN | 4.4.0 |
| `fansi` | 1.0.6 | CRAN | 4.4.0 |
| `farver` | 2.1.2 | CRAN | 4.4.0 |
| `fastcluster` | 1.3.0 | CRAN | 4.4.1 |
| `fastmap` | 1.2.0 | CRAN | 4.4.0 |
| `fastmatch` | 1.1.4 | CRAN | 4.4.0 |
| `fgsea` | 1.30.0 | Bioconductor | 4.4.0 |
| `filelock` | 1.0.3 | CRAN | 4.4.0 |
| `flowCore` | 2.16.0 | Bioconductor | 4.4.0 |
| `fontawesome` | 0.5.3 | CRAN | 4.4.1 |
| `forcats` | 1.0.1 | CRAN | 4.4.1 |
| `foreach` | 1.5.2 | CRAN | 4.4.0 |
| `forecast` | 9.0.2 | CRAN | 4.4.3 |
| `foreign` | 0.8.86 | CRAN | 4.4.1 |
| `formatR` | 1.14 | CRAN | 4.4.0 |
| `Formula` | 1.2.5 | CRAN | 4.4.0 |
| `fracdiff` | 1.5.3 | CRAN | 4.4.1 |
| `fs` | 1.6.7 | CRAN | 4.4.3 |
| `futile.logger` | 1.4.9 | CRAN | 4.4.3 |
| `futile.options` | 1.0.1 | CRAN | 4.4.0 |
| `gargle` | 1.5.2 | CRAN | 4.4.0 |
| `gdsfmt` | 1.42.1 | Bioconductor 3.20 | 4.4.2 |
| `generics` | 0.1.4 | CRAN | 4.4.1 |
| `GenomeInfoDb` | 1.40.1 | Bioconductor | 4.4.0 |
| `GenomeInfoDbData` | 1.2.12 | Bioconductor | 4.4.1 |
| `GenomicAlignments` | 1.40.0 | Bioconductor | 4.4.0 |
| `GenomicFeatures` | 1.56.0 | Bioconductor | 4.4.0 |
| `GenomicRanges` | 1.58.0 | Bioconductor 3.20 | 4.4.1 |
| `geomtextpath` | 0.2.0 | CRAN | 4.4.1 |
| `gert` | 2.1.4 | CRAN | 4.4.1 |
| `GetoptLong` | 1.1.0 | CRAN | 4.4.3 |
| `ggbeeswarm` | 0.7.2 | CRAN | 4.4.0 |
| `ggfittext` | 0.10.3 | CRAN | 4.4.3 |
| `ggforce` | 0.5.0 | CRAN | 4.4.1 |
| `ggfun` | 0.1.7 | CRAN | 4.4.1 |
| `ggh4x` | 0.3.1 | CRAN | 4.4.1 |
| `ggnewscale` | 0.5.0 | CRAN | 4.4.0 |
| `ggpattern` | 1.3.1 | CRAN | 4.4.3 |
| `ggplot2` | 4.0.2 | CRAN | 4.4.3 |
| `ggplotify` | 0.1.2 | CRAN | 4.4.0 |
| `ggpubr` | 0.6.3 | CRAN | 4.4.3 |
| `ggraph` | 2.2.2 | CRAN | 4.4.1 |
| `ggrastr` | 1.0.2 | CRAN | 4.4.0 |
| `ggrepel` | 0.9.8 | CRAN | 4.4.3 |
| `ggsankey` | 0.0.99999 | github:davidsjoberg/ggsankey@b675d0d | 4.4.1 |
| `ggsci` | 4.2.0 | CRAN | 4.4.3 |
| `ggsignif` | 0.6.4 | CRAN | 4.4.0 |
| `ggtext` | 0.1.2 | CRAN | 4.4.0 |
| `ggtree` | 3.12.0 | Bioconductor | 4.4.0 |
| `ggVennDiagram` | 1.5.2 | CRAN | 4.4.0 |
| `gh` | 1.4.1 | CRAN | 4.4.0 |
| `gitcreds` | 0.1.2 | CRAN | 4.4.0 |
| `glmnet` | 4.1.8 | CRAN | 4.4.0 |
| `GlobalOptions` | 0.1.3 | CRAN | 4.4.3 |
| `glue` | 1.8.1 | CRAN | 4.4.1 |
| `GO.db` | 3.19.1 | Bioconductor | 4.4.1 |
| `googledrive` | 2.1.1 | CRAN | 4.4.0 |
| `googlesheets4` | 1.1.1 | CRAN | 4.4.0 |
| `GOSemSim` | 2.30.2 | Bioconductor | 4.4.1 |
| `GPArotation` | 2024.3.1 | CRAN | 4.4.0 |
| `gplots` | 3.3.0 | CRAN | 4.4.3 |
| `gprofiler2` | 0.2.3 | CRAN | 4.4.0 |
| `graph` | 1.82.0 | Bioconductor | 4.4.0 |
| `graphics` | 4.4.1 | base | 4.4.1 |
| `graphlayouts` | 1.2.2 | CRAN | 4.4.1 |
| `grDevices` | 4.4.1 | base | 4.4.1 |
| `grid` | 4.4.1 | base | 4.4.1 |
| `gridExtra` | 2.3 | CRAN | 4.4.1 |
| `gridGraphics` | 0.5.1 | CRAN | 4.4.0 |
| `gridpattern` | 1.3.1 | CRAN | 4.4.1 |
| `gridtext` | 0.1.5 | CRAN | 4.4.0 |
| `GSEABase` | 1.66.0 | Bioconductor | 4.4.0 |
| `gson` | 0.1.0 | CRAN | 4.4.0 |
| `gtable` | 0.3.6 | CRAN | 4.4.1 |
| `gtools` | 3.9.5 | CRAN | 4.4.0 |
| `haven` | 2.5.4 | CRAN | 4.4.0 |
| `here` | 1.0.2 | CRAN | 4.4.1 |
| `highr` | 0.12 | CRAN | 4.4.3 |
| `Hmisc` | 5.2.0 | CRAN | 4.4.1 |
| `hms` | 1.1.4 | CRAN | 4.4.1 |
| `htmlTable` | 2.4.3 | CRAN | 4.4.0 |
| `htmltools` | 0.5.9 | CRAN | 4.4.3 |
| `htmlwidgets` | 1.6.4 | CRAN | 4.4.0 |
| `httpuv` | 1.6.15 | CRAN | 4.4.0 |
| `httr` | 1.4.8 | CRAN | 4.4.3 |
| `httr2` | 1.0.7 | CRAN | 4.4.1 |
| `ids` | 1.0.1 | CRAN | 4.4.0 |
| `igraph` | 2.2.1 | CRAN | 4.4.1 |
| `impute` | 1.80.0 | Bioconductor 3.20 | 4.4.1 |
| `ini` | 0.3.1 | CRAN | 4.4.0 |
| `inspectdf` | 0.0.12.1 | CRAN | 4.4.1 |
| `IRanges` | 2.40.1 | Bioconductor 3.20 | 4.4.2 |
| `isoband` | 0.3.0 | CRAN | 4.4.3 |
| `iterators` | 1.0.14 | CRAN | 4.4.0 |
| `jomo` | 2.7.6 | CRAN | 4.4.0 |
| `jpeg` | 0.1.11 | CRAN | 4.4.1 |
| `jquerylib` | 0.1.4 | CRAN | 4.4.0 |
| `jsonlite` | 2.0.0 | CRAN | 4.4.1 |
| `KEGGREST` | 1.44.1 | Bioconductor | 4.4.0 |
| `KernSmooth` | 2.23.24 | CRAN | 4.4.1 |
| `knitr` | 1.51 | CRAN | 4.4.3 |
| `labeling` | 0.4.3 | CRAN | 4.4.0 |
| `lambda.r` | 1.2.4 | CRAN | 4.4.0 |
| `later` | 1.4.1 | CRAN | 4.4.1 |
| `latex2exp` | 0.9.8 | CRAN | 4.4.3 |
| `lattice` | 0.22.6 | CRAN | 4.4.1 |
| `lazyeval` | 0.2.2 | CRAN | 4.4.0 |
| `lifecycle` | 1.0.5 | CRAN | 4.4.3 |
| `limma` | 3.62.2 | Bioconductor | 4.4.2 |
| `litedown` | 0.9 | CRAN | 4.4.3 |
| `lme4` | 1.1.38 | CRAN | 4.4.1 |
| `lmerTest` | 3.2.1 | CRAN | 4.4.3 |
| `lmtest` | 0.9.40 | CRAN | 4.4.0 |
| `locfit` | 1.5.9.12 | CRAN | 4.4.1 |
| `lubridate` | 1.9.5 | CRAN | 4.4.3 |
| `magrittr` | 2.0.5 | CRAN | 4.4.1 |
| `markdown` | 2.0 | CRAN | 4.4.1 |
| `MASS` | 7.3.60.2 | CRAN | 4.4.1 |
| `mathjaxr` | 1.6.0 | CRAN | 4.4.0 |
| `Matrix` | 1.7.0 | CRAN | 4.4.1 |
| `MatrixGenerics` | 1.16.0 | Bioconductor | 4.4.0 |
| `MatrixModels` | 0.5.4 | CRAN | 4.4.1 |
| `matrixStats` | 1.5.0 | CRAN | 4.4.1 |
| `memoise` | 2.0.1 | CRAN | 4.4.0 |
| `metap` | 1.11 | CRAN | 4.4.0 |
| `methods` | 4.4.1 | base | 4.4.1 |
| `Mfuzz` | 2.72.0 | Bioconductor | 4.4.1 |
| `mgcv` | 1.9.1 | CRAN | 4.4.1 |
| `mice` | 3.17.0 | CRAN | 4.4.1 |
| `microbenchmark` | 1.5.0 | CRAN | 4.4.1 |
| `mime` | 0.13 | CRAN | 4.4.1 |
| `miniUI` | 0.1.1.1 | CRAN | 4.4.0 |
| `minqa` | 1.2.8 | CRAN | 4.4.0 |
| `mitml` | 0.4.5 | CRAN | 4.4.0 |
| `mnormt` | 2.1.1 | CRAN | 4.4.0 |
| `modelr` | 0.1.11 | CRAN | 4.4.0 |
| `MotrpacBicQC` | 1.7.0 | github:MoTrPAC/MotrpacBicQC@d873a93 | 4.4.1 |
| `MotrpacHumanPreSuspensionAnalysis` | 0.2.4 | local source install | 4.4.1 |
| `MotrpacHumanPreSuspensionData` | 0.0.1.102 | local source install | 4.4.1 |
| `MotrpacRatTraining6mo` | 1.6.5 | github:MoTrPAC/MotrpacRatTraining6mo@0f38357 | 4.4.1 |
| `MotrpacRatTraining6moData` | 2.0.0 | github:MoTrPAC/MotrpacRatTraining6moData@874814b | 4.4.1 |
| `msigdbr` | 25.1.1 | CRAN | 4.4.1 |
| `multcomp` | 1.4.26 | CRAN | 4.4.0 |
| `multtest` | 2.60.0 | Bioconductor | 4.4.0 |
| `mutoss` | 0.1.13 | CRAN | 4.4.0 |
| `mvtnorm` | 1.3.6 | CRAN | 4.4.3 |
| `naniar` | 1.1.0 | CRAN | 4.4.0 |
| `nlme` | 3.1.164 | CRAN | 4.4.1 |
| `nloptr` | 2.2.1 | CRAN | 4.4.1 |
| `nnet` | 7.3.19 | CRAN | 4.4.1 |
| `norm` | 1.0.11.1 | CRAN | 4.4.0 |
| `nortest` | 1.0.4 | CRAN | 4.4.0 |
| `numDeriv` | 2016.8.1.1 | CRAN | 4.4.0 |
| `openssl` | 2.3.5 | CRAN | 4.4.3 |
| `openxlsx` | 4.2.8 | CRAN | 4.4.1 |
| `ordinal` | 2023.12.4.1 | CRAN | 4.4.1 |
| `org.Hs.eg.db` | 3.19.1 | Bioconductor | 4.4.1 |
| `pak` | 0.8.0 | CRAN | 4.4.1 |
| `pan` | 1.9 | CRAN | 4.4.0 |
| `parallel` | 4.4.1 | base | 4.4.1 |
| `patchwork` | 1.3.2 | CRAN | 4.4.1 |
| `pbkrtest` | 0.5.5 | CRAN | 4.4.1 |
| `pheatmap` | 1.0.12 | CRAN | 4.4.0 |
| `pillar` | 1.11.1 | CRAN | 4.4.1 |
| `pkgbuild` | 1.4.8 | CRAN | 4.4.1 |
| `pkgconfig` | 2.0.3 | CRAN | 4.4.0 |
| `pkgdown` | 2.1.1 | CRAN | 4.4.1 |
| `pkgload` | 1.5.2 | CRAN | 4.4.1 |
| `PLIER` | 0.99.0 | github:wgmao/PLIER@fe4e9b2 | 4.4.1 |
| `plotly` | 4.11.0 | CRAN | 4.4.1 |
| `plotrix` | 3.8.6 | github:plotrix/plotrix@0d4c2b0 | 4.4.1 |
| `plyr` | 1.8.9 | CRAN | 4.4.0 |
| `png` | 0.1.9 | CRAN | 4.4.3 |
| `polyclip` | 1.10.7 | CRAN | 4.4.0 |
| `polynom` | 1.4.1 | CRAN | 4.4.0 |
| `praise` | 1.0.0 | CRAN | 4.4.0 |
| `preprocessCore` | 1.68.0 | Bioconductor 3.20 | 4.4.1 |
| `prettyunits` | 1.2.0 | CRAN | 4.4.0 |
| `processx` | 3.8.6 | CRAN | 4.4.1 |
| `profvis` | 0.4.0 | CRAN | 4.4.1 |
| `progress` | 1.2.3 | CRAN | 4.4.0 |
| `promises` | 1.3.2 | CRAN | 4.4.1 |
| `ProtGenerics` | 1.38.0 | Bioconductor 3.20 | 4.4.1 |
| `proxy` | 0.4.27 | CRAN | 4.4.0 |
| `ps` | 1.9.1 | CRAN | 4.4.1 |
| `psych` | 2.4.6.26 | CRAN | 4.4.0 |
| `purrr` | 1.2.1 | CRAN | 4.4.3 |
| `qqconf` | 1.3.2 | CRAN | 4.4.0 |
| `quantreg` | 6.1 | CRAN | 4.4.1 |
| `qvalue` | 2.36.0 | Bioconductor | 4.4.0 |
| `R.methodsS3` | 1.8.2 | CRAN | 4.4.0 |
| `R.oo` | 1.27.0 | CRAN | 4.4.1 |
| `R.utils` | 2.12.3 | CRAN | 4.4.0 |
| `r2r` | 0.1.2 | CRAN | 4.4.1 |
| `R6` | 2.6.1 | CRAN | 4.4.1 |
| `ragg` | 1.3.3 | CRAN | 4.4.1 |
| `randomForest` | 4.7.1.2 | CRAN | 4.4.1 |
| `rappdirs` | 0.3.4 | CRAN | 4.4.3 |
| `rbibutils` | 2.4.1 | CRAN | 4.4.3 |
| `rcmdcheck` | 1.4.0 | CRAN | 4.4.0 |
| `RColorBrewer` | 1.1.3 | CRAN | 4.4.0 |
| `Rcpp` | 1.1.1 | CRAN | 4.4.3 |
| `RcppArmadillo` | 15.2.4.1 | CRAN | 4.4.3 |
| `RcppEigen` | 0.3.4.0.2 | CRAN | 4.4.1 |
| `RcppTOML` | 0.2.3 | CRAN | 4.4.1 |
| `RCurl` | 1.98.1.17 | CRAN | 4.4.1 |
| `Rdpack` | 2.6.6 | CRAN | 4.4.3 |
| `readr` | 2.2.0 | CRAN | 4.4.3 |
| `readxl` | 1.4.3 | CRAN | 4.4.0 |
| `reformulas` | 0.4.4 | CRAN | 4.4.3 |
| `regioneR` | 1.38.0 | Bioconductor 3.20 | 4.4.1 |
| `remaCor` | 0.0.20 | CRAN | 4.4.1 |
| `rematch` | 2.0.0 | CRAN | 4.4.0 |
| `rematch2` | 2.1.2 | CRAN | 4.4.0 |
| `remotes` | 2.5.0 | CRAN | 4.4.0 |
| `renv` | 1.0.11 | CRAN | 4.4.1 |
| `reprex` | 2.1.1 | CRAN | 4.4.0 |
| `reshape` | 0.8.10 | CRAN | 4.4.1 |
| `reshape2` | 1.4.5 | CRAN | 4.4.1 |
| `restfulr` | 0.0.15 | CRAN | 4.4.0 |
| `reticulate` | 1.43.0 | CRAN | 4.4.1 |
| `rhdf5` | 2.48.0 | Bioconductor | 4.4.0 |
| `rhdf5filters` | 1.16.0 | Bioconductor | 4.4.0 |
| `Rhdf5lib` | 1.26.0 | Bioconductor | 4.4.0 |
| `RhpcBLASctl` | 0.23.42 | CRAN | 4.4.0 |
| `Rhtslib` | 3.0.0 | Bioconductor | 4.4.0 |
| `rjson` | 0.2.23 | CRAN | 4.4.1 |
| `rlang` | 1.3.0 | CRAN | 4.4.1 |
| `rmarkdown` | 2.30 | CRAN | 4.4.1 |
| `Rmisc` | 1.5.1 | CRAN | 4.4.0 |
| `robustbase` | 0.99.4.1 | CRAN | 4.4.1 |
| `roxygen2` | 8.0.0 | CRAN | 4.4.1 |
| `rpart` | 4.1.23 | CRAN | 4.4.1 |
| `rprojroot` | 2.1.1 | CRAN | 4.4.1 |
| `RProtoBufLib` | 2.16.0 | Bioconductor | 4.4.0 |
| `Rsamtools` | 2.20.0 | Bioconductor | 4.4.0 |
| `RSpectra` | 0.16.2 | CRAN | 4.4.0 |
| `RSQLite` | 2.4.6 | CRAN | 4.4.3 |
| `rstatix` | 0.7.3 | CRAN | 4.4.1 |
| `rstudioapi` | 0.17.1 | CRAN | 4.4.1 |
| `rsvd` | 1.0.5 | CRAN | 4.4.0 |
| `rtracklayer` | 1.64.0 | Bioconductor | 4.4.0 |
| `Rtsne` | 0.17 | CRAN | 4.4.0 |
| `rversions` | 2.1.2 | CRAN | 4.4.0 |
| `rvest` | 1.0.4 | CRAN | 4.4.0 |
| `s2` | 1.1.9 | CRAN | 4.4.1 |
| `S4Arrays` | 1.6.0 | Bioconductor 3.20 | 4.4.1 |
| `S4Vectors` | 0.44.0 | Bioconductor 3.20 | 4.4.1 |
| `S7` | 0.2.1 | CRAN | 4.4.3 |
| `sandwich` | 3.1.1 | CRAN | 4.4.1 |
| `sass` | 0.4.10 | CRAN | 4.4.1 |
| `scales` | 1.4.0 | CRAN | 4.4.1 |
| `scatterpie` | 0.2.4 | CRAN | 4.4.1 |
| `selectr` | 0.4.2 | CRAN | 4.4.0 |
| `sessioninfo` | 1.2.2 | CRAN | 4.4.0 |
| `sf` | 1.0.21 | CRAN | 4.4.1 |
| `shades` | 1.4.0 | CRAN | 4.4.0 |
| `shadowtext` | 0.1.4 | CRAN | 4.4.0 |
| `shape` | 1.4.6.1 | CRAN | 4.4.0 |
| `shiny` | 1.9.1 | CRAN | 4.4.0 |
| `sn` | 2.1.1 | CRAN | 4.4.0 |
| `snow` | 0.4.4 | CRAN | 4.4.0 |
| `SNPRelate` | 1.40.0 | Bioconductor 3.20 | 4.4.1 |
| `sourcetools` | 0.1.7.1 | CRAN | 4.4.0 |
| `SparseArray` | 1.6.2 | Bioconductor 3.20 | 4.4.2 |
| `SparseM` | 1.84.2 | CRAN | 4.4.0 |
| `splines` | 4.4.1 | base | 4.4.1 |
| `statmod` | 1.5.1 | CRAN | 4.4.1 |
| `stats` | 4.4.1 | base | 4.4.1 |
| `stats4` | 4.4.1 | base | 4.4.1 |
| `stringi` | 1.8.7 | CRAN | 4.4.1 |
| `stringr` | 1.6.0 | CRAN | 4.4.1 |
| `SummarizedExperiment` | 1.34.0 | Bioconductor | 4.4.0 |
| `survival` | 3.6.4 | CRAN | 4.4.1 |
| `sys` | 3.4.3 | CRAN | 4.4.1 |
| `systemfonts` | 1.3.1 | CRAN | 4.4.1 |
| `table.glue` | 0.0.5 | CRAN | 4.4.1 |
| `tcltk` | 4.4.1 | base | 4.4.1 |
| `testthat` | 3.3.2 | CRAN | 4.4.3 |
| `textshaping` | 0.4.0 | CRAN | 4.4.0 |
| `TFisher` | 0.2.0 | CRAN | 4.4.0 |
| `TH.data` | 1.1.2 | CRAN | 4.4.0 |
| `tibble` | 3.3.1 | CRAN | 4.4.3 |
| `tidygraph` | 1.3.1 | CRAN | 4.4.0 |
| `tidyr` | 1.3.2 | CRAN | 4.4.3 |
| `tidyselect` | 1.2.1 | CRAN | 4.4.0 |
| `tidytree` | 0.4.6 | CRAN | 4.4.0 |
| `tidyverse` | 2.0.0 | CRAN | 4.4.0 |
| `timechange` | 0.4.0 | CRAN | 4.4.3 |
| `timeDate` | 4052.112 | CRAN | 4.4.3 |
| `tinytex` | 0.58 | CRAN | 4.4.3 |
| `tkWidgets` | 1.82.0 | Bioconductor | 4.4.0 |
| `TMSig` | 1.6.0 | Bioconductor | 4.4.1 |
| `tools` | 4.4.1 | base | 4.4.1 |
| `treeio` | 1.28.0 | Bioconductor | 4.4.0 |
| `tweenr` | 2.0.3 | CRAN | 4.4.0 |
| `TxDb.Hsapiens.UCSC.hg19.knownGene` | 3.2.2 | Bioconductor 3.20 | 4.4.1 |
| `txdbmaker` | 1.2.1 | Bioconductor 3.20 | 4.4.2 |
| `tzdb` | 0.5.0 | CRAN | 4.4.1 |
| `ucminf` | 1.2.2 | CRAN | 4.4.0 |
| `UCSC.utils` | 1.0.0 | Bioconductor | 4.4.0 |
| `umap` | 0.2.10.0 | CRAN | 4.4.0 |
| `units` | 0.8.7 | CRAN | 4.4.1 |
| `UpSetR` | 1.4.0 | CRAN | 4.4.0 |
| `urca` | 1.3.4 | CRAN | 4.4.0 |
| `urlchecker` | 1.0.1 | CRAN | 4.4.0 |
| `usethis` | 3.1.0 | CRAN | 4.4.1 |
| `utf8` | 1.2.6 | CRAN | 4.4.1 |
| `utils` | 4.4.1 | base | 4.4.1 |
| `uuid` | 1.2.1 | CRAN | 4.4.0 |
| `variancePartition` | 1.38.1 | Bioconductor | 4.4.1 |
| `vctrs` | 0.7.3 | CRAN | 4.4.1 |
| `venn` | 1.12 | CRAN | 4.4.1 |
| `VennDiagram` | 1.7.3 | CRAN | 4.4.0 |
| `vipor` | 0.4.7 | CRAN | 4.4.1 |
| `viridis` | 0.6.5 | CRAN | 4.4.0 |
| `viridisLite` | 0.4.3 | CRAN | 4.4.3 |
| `visdat` | 0.6.0 | CRAN | 4.4.0 |
| `visNetwork` | 2.1.2 | CRAN | 4.4.0 |
| `vroom` | 1.7.0 | CRAN | 4.4.3 |
| `waldo` | 0.6.2 | CRAN | 4.4.1 |
| `WGCNA` | 1.73 | CRAN | 4.4.1 |
| `whisker` | 0.4.1 | CRAN | 4.4.0 |
| `widgetTools` | 1.82.0 | Bioconductor | 4.4.0 |
| `withr` | 3.0.3 | CRAN | 4.4.1 |
| `wk` | 0.9.4 | CRAN | 4.4.1 |
| `writexl` | 1.5.4 | CRAN | 4.4.1 |
| `xfun` | 0.57 | CRAN | 4.4.3 |
| `XML` | 3.99.0.22 | CRAN | 4.4.3 |
| `xml2` | 1.5.2 | CRAN | 4.4.3 |
| `xopen` | 1.0.1 | CRAN | 4.4.0 |
| `xtable` | 1.8.8 | CRAN | 4.4.3 |
| `XVector` | 0.44.0 | Bioconductor | 4.4.0 |
| `yaml` | 2.3.12 | CRAN | 4.4.3 |
| `yulab.utils` | 0.2.0 | CRAN | 4.4.1 |
| `zip` | 2.3.3 | CRAN | 4.4.1 |
| `zlibbioc` | 1.50.0 | Bioconductor | 4.4.0 |
| `zoo` | 1.8.15 | CRAN | 4.4.3 |

</details>

## R session

```
R version:       R version 4.4.1 (2024-06-14)
Platform:        aarch64-apple-darwin20
Bioconductor:    3.19
Packages listed: 436 (118 directly declared, 318 transitive)
Referenced in source but not installed (may include false positives where a variable was scanned as a package name): BSgenome.Hsapiens.UCSC.hg38, EnsDb.Hsapiens.v86, FactoMineR, ggnetwork, magical, motifmatchr, MotrpacHumanPreSuspension, org.Mm.eg.db, Orthology.eg.db, PMA, TFBSTools

Declared dependencies per repo:
  MotrpacHumanPreSuspensionData: 14 (DESCRIPTION)
  MotrpacHumanPreSuspensionAnalysis: 31 (DESCRIPTION)
  MotrpacPreSuspensionAcute: 113 (source scan)
  motrpac-human-presuspension-repro: 36 (source scan)

.libPaths():
  /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library

sessionInfo():
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS 26.5.2

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] BiocManager_1.30.26 compiler_4.4.1      tools_4.4.1        
```

---
Captured by `config/capture_environment.sh` on 2026-07-30 09:46:30 CDT (-0500).
