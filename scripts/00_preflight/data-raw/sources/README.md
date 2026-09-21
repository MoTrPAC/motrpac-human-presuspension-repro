# Non-GMT builder raw sources

Raw inputs for the non-GMT leaf builders in `../`. Provenance below is taken
from the original builder-script comments; where a file is a working-repo asset
rather than an external download, no download date exists.

| Source file | Provenance | Downloaded / version | Builder | License |
|-------------|-----------|----------------------|---------|---------|
| `covariates_pre_cawg.csv` | precovid working repo `library/covariates_pre_cawg.csv` (managed by Christopher Jin) | repo asset — no download date | `../COVARIATES_FILE.R` | internal |

Note: `covariates_pre_cawg.csv` is not an external download; upstream
`COVARIATES_FILE.R` reads it from `precovid_repo_path` in `~/config.json`. It is
vendored here so the builder runs without that config entry.

## `ensembl_v105/` — local Ensembl release 105 annotation cache
Pure external reference data (Ensembl v105 TxDb + gene-attribute table) with no
dependency on any study/quant-id data, so it lives at the preflight level. Built
once by `ensembl_v105/build_ensembl_v105_cache.R` and consumed OFFLINE by the
Stage 1 QC-norm annotation stems (via `lib/qc_helpers.R`), replacing the flaky
live biomaRt/txdbmaker calls. See `ensembl_v105/README.md` for full provenance.

## `refmet/` — local RefMet / KEGG annotation snapshot
Frozen responses from the Metabolomics Workbench RefMet batch endpoint and
`KEGGREST::keggList("compound")`, built once by `refmet/build_refmet_cache.R` and
consumed OFFLINE by the Stage 1 metabolomics QC-norm stem (via `lib/qc_helpers.R`),
replacing the two live queries that made `refmet_name`, `refmet_id` and `kegg_id`
depend on the run date. RefMet is a rolling database with no releases, so the
snapshot's build date is its version. Unlike `ensembl_v105/`, the snapshot is keyed
by the names the study asks about — the endpoint resolves lab-supplied synonyms
server-side, so a full-database dump could not stand in for it. Refreshing therefore
needs the (gated) raw metabolomics metadata to re-derive that name universe, which
is why this folder is only *mostly* study-independent. See `refmet/README.md`.

---

# Gated builder raw sources

Raw inputs for the gated leaf builders in `../`. These are **consortium-gated
MoTrPAC data** (sample/vial-level) — do not redistribute outside an approved
data-access agreement.

## `removed_samples/` — for `../OUTLIERS.R`
21 `*_removed-samples_*.txt` files downloaded **2026-07-27** from the production
data hub:
`gs://motrpac-data-hub/analysis/human-precovid-sed-adu/v1.3/{omics_group}/metadata/`.
Each filename encodes the tissue code (`t0x-…`) and ome, which `OUTLIERS.R` reads
via the vendored `.find_tissue` / `.find_ome`. Column schemas vary across files
(`vialLabel`/`reason`, `sample`/`Tissue`/`reason`, `Dataset`/`tissue`/`ome`/`Sample`/`PC`);
the builder normalizes them (`sample`/`Sample`→`vialLabel`, `PC`→`reason`).
To refresh: re-run the `gsutil cp` of `**removed-samples*.txt` from the v1.3 bucket.

## `differential_splicing.rds` — for `../SPLICING_DA.R`  (NOT vendored)
Differential alternative-splicing results (FDR < 0.05), supplied by **Zidong
Zhang, MSSM (zidong.zhang@mssm.edu)** — not downloadable from a bucket. Place the
file here by hand to build `SPLICING_DA`; absent that, the builder errors with
instructions.
