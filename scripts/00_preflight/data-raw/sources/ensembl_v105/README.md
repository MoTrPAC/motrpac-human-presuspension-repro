# Local Ensembl v105 annotation cache

Local, offline snapshot of **Ensembl release 105** used by the Stage 1 QC-norm
stems to annotate features. The live Ensembl/biomaRt connection is flaky
(intermittent `curl_fetch_memory` errors, `Fetch exons and CDS ... Error: [0]`,
"try a mirror site"), which repeatedly broke annotation-dependent stems (ATAC,
methylcap, prot-ol/pr/ph, clinical). Caching v105 once makes every run offline
and deterministic; the only script that contacts Ensembl is the builder below.

## Artifacts

| File | Provenance | Consumed by |
|------|-----------|-------------|
| `txdb_hsapiens_ensembl_v105.sqlite` | `txdbmaker::makeTxDbFromEnsembl("Homo Sapiens", release=105)` → `AnnotationDbi::saveDb` | ChIPseeker peak annotation in `pre_cawg_get_peak_annotations_hs()` (atac, methylcap) |
| `ensembl_v105_gene_attributes.rds` | one `biomaRt::getBM` (release 105, all genes) of `ensembl_gene_id, entrezgene_id, external_gene_name, uniprotswissprot` | every `getBM` in the stems, via `ensembl_v105_getBM()` in `lib/qc_helpers.R` |

The gene-attribute table is the **union** of every stem's requested attributes;
each stem filters it by `ensembl_gene_id`, `external_gene_name`, or
`uniprotswissprot` locally (no network).

## Build / refresh

`txdb_hsapiens_ensembl_v105.sqlite` is gitignored (175 MB, over GitHub's file
limit), so it is distributed through the bucket instead: step 06's
`stage_ensembl_txdb` stem copies it into `staging/freeze/resources/`, Stage 2
uploads it, and a fresh clone gets it back with

```bash
make sources                      # downloads it from resources/
```

The published copy carries a `_v<version>` suffix from v2.0 on (its entry is under
`resources` in `config/file_versions.json`); this local cache keeps the bare name,
since it is an input to the pipeline rather than a release artifact.
`copy_from_source.sh` globs both spellings and takes the highest, so a bucket still
holding the v1.4-era unversioned copy fetches just the same.

Rebuilding from Ensembl is the fallback for when `resources/` has no copy yet:

```bash
Rscript build_ensembl_v105_cache.R   # FORCE=1 to rebuild artifacts that exist
```

Each artifact is skipped when it already exists, so re-running the builder is cheap.
`ensembl_v105_gene_attributes.rds` is small enough to be committed, so it is left
alone unless `FORCE=1`.

Requires network access to Ensembl (release 105) and the `txdbmaker`, `biomaRt`,
`AnnotationDbi` packages. The builder retries the (flaky) fetches a bounded number
of times and errors loudly if Ensembl stays unreachable — re-run when it is up.
**Built: 2026-07-27** (Ensembl release 105).

## Reproducibility note

The three epigen/transcriptomics stems (atac, methylcap, transcriptomics) were
already pinned to v105. The four proteomics/clinical stems (prot-ol, prot-pr,
prot-ph, clinical) previously used `biomaRt::useMart(...)` **unpinned** (latest
Ensembl); switching them to this v105 cache pins them to v105 — more reproducible,
but their feature-metadata outputs were re-diffed against `staging_20260720` to
confirm they did not regress.
