# Preflight data objects (`data-raw/`)

Self-contained builders for the **leaf data objects** — the ones the inventory
(`docs/data_objects.tsv`) marks `deps=-`, i.e. built only from external
assets with no upstream pipeline dependency. They are vendored + adapted from
`MotrpacHumanPreSuspensionAnalysis/data-raw/` so this repo can regenerate them
without installing the Analysis package. Outputs are split R-package style:

- `../data/`   — `.rds` data objects **only** (COLORS, COVARIATES_FILE,
  OME_TISSUE_CODE, OUTLIERS, SPLICING_DA).
- `../output/` — everything else: the `*.gmt.gz` GMTs and
  `quantid_bucket_files.json` (the quant-id catalog, which `OME_TISSUE_CODE.R`
  reads as input).

Run all of them:

```bash
bash scripts/00_preflight/data-raw/build_data_objects.sh      # or: make data-objects
bash scripts/00_preflight/data-raw/build_data_objects.sh GMT_ # subset by name match
```

## What gets built

All builders live directly in `data-raw/` (no gated/non-gated split).

| Object(s) | Builder | Source |
|-----------|---------|--------|
| `HUMAN_*` palettes (6 rds) | `COLORS_ABBREVIATIONS.R` | hardcoded literals |
| `COVARIATES_FILE.rds` | `COVARIATES_FILE.R` | `sources/covariates_pre_cawg.csv` |
| `QUANTID_BUCKET_FILES.rds` | `QUANTID_BUCKET_FILES.R` | `../output/quantid_bucket_files.json` |
| `OME_TISSUE_CODE.rds` | `OME_TISSUE_CODE.R` | **`../output/quantid_bucket_files.json`** (not a live `gsutil ls`) |
| `cellmarker.v2024…gmt.gz` | `gmt_processing/GMT_CellMarker.R` | `gmt_processing/sources/CellMarker_2024.txt.gz` |
| `mitocarta3.0…gmt.gz` | `gmt_processing/GMT_MitoCarta.R` | `…/Human.MitoCarta3.0.xls` |
| `phosphositeplus…gmt.gz` | `gmt_processing/GMT_PSP_kinase.R` | `…/Kinase_Substrate_Dataset.gz` |
| `ptmsigdb.v2.0…gmt.gz` | `gmt_processing/GMT_PTMSigDB.R` | `…/ptm.sig.db…v2.0.0.gmt.gz` |
| `…refmet.2024.08.07…gmt.gz` | `gmt_processing/GMT_RefMet_subclass.R` | frozen snapshot in `gmt_processing/prebuilt/`* |
| `c2/c5…v2023.2…gmt.gz` | (copied by the driver) | `gmt_processing/prebuilt/` (MSigDB, login-gated download)* |
| `OUTLIERS.rds` | `OUTLIERS.R` | `sources/removed_samples/` (23 files) |
| `SPLICING_DA.rds` | `SPLICING_DA.R` | `sources/differential_splicing.rds` (from MSSM) |

\* Not rebuilt locally: MSigDB C2/C5 are login-gated downloads and RefMet needs
`HUMAN_FEATURE_TO_GENE` (upstream, non-leaf) plus a manual web-form submission;
their committed artifacts are vendored under `gmt_processing/prebuilt/`. The
upstream RefMet build is kept, commented out, in `GMT_RefMet_subclass.R`.

`OUTLIERS.R` reads removed-samples files vendored under `sources/removed_samples/`
(rather than fetching live), and `SPLICING_DA` needs `sources/differential_splicing.rds`
from MSSM (SKIPs if absent). Both underlying files are **consortium data**
(sample/vial-level) — don't redistribute; see `sources/README.md`.

## Provenance & dates
Download dates / versions for every raw input are recorded in the per-folder
READMEs: `gmt_processing/sources/`, `gmt_processing/prebuilt/`, and `sources/`.
Dates come from the original builder-script comments.

## Fidelity
The rebuilt GMTs were checked against the committed
`MotrpacHumanPreSuspensionAnalysis/.../gmt_files/`: MitoCarta is byte-identical
(decompressed), and CellMarker / PSP / PTMSigDB have **identical set membership**
(only within-line member ordering differs, which is irrelevant for gene sets).

## Keep in sync
These are **vendored copies**. If an upstream builder in the Analysis package
changes, update the copy here. Known upstream bugs fixed in these copies: the
`gzip(file)` / `.writeGMT(path=file)` undefined-variable bugs in
`GMT_MitoCarta.R`, `GMT_PSP_kinase.R`, `GMT_PTMSigDB.R`, `GMT_RefMet_subclass.R`.
