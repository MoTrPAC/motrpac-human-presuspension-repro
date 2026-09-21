# Step 05 (stage_metabolomics_cvs) raw sources

Same convention as `../../06_generate_qc_norm/sources/` and `../../09_build_da/sources/`:
inputs this pipeline does not regenerate, vendored so the step runs without reaching
outside the repo.

## `METABOLOMICS_CVS.rda`

The per-metabolite coefficient-of-variation table, one row per
tissue × assay × site × feature (4000 rows, 9 columns).

### Origin

The CVs are calculated by

```
precovid-analyses/QC/human-precovid-sed-adu_all_metabolomics_qc-cvs.Rmd
```

("Metabolomics Human-precovid: CV calculations", David Jimenez-Morales and Gayatri
Iyer), which builds a `metab_cvs` table. Every copy of this object downstream is a
version of that single output:

| Artifact | How it derives |
|----------|----------------|
| `.../c1.3/resources/motrpac_human-precovid_metabolite-cv_v1.2.txt` | the Rmd's table, published to the production bucket |
| `MotrpacHumanPreSuspensionAnalysis/data/METABOLOMICS_CVS.rda` | the Analysis package's `data-raw/METABOLOMICS_CVS.R` downloads that bucket `.txt` and `use_data()`s it — it computes no CVs of its own |
| `sources/METABOLOMICS_CVS.rda` (this file) | copied from that package `data/` |
| `staging/freeze/resources/motrpac_human-precovid_metabolite-cv_v1.2.txt` | `../METABOLOMICS_CVS_resource.R` re-emits this file |

**Provenance:** copied **2026-07-29** from
`MotrpacHumanPreSuspensionAnalysis/data/METABOLOMICS_CVS.rda`
(sha1 `f876cb1e173ca97c700003f2978827dab06e8d00`).

Vendored rather than re-copied from that checkout on every run, so step 05 does not
depend on the Analysis package being present or on whatever state its `data/` is in.
`build.sh` copies this into `scripts/10_build_data/data/`, and
`../METABOLOMICS_CVS_resource.R` also emits it as
`staging/freeze/resources/motrpac_human-precovid_metabolite-cv_v1.2.txt`.

**This is a stage, not a local build.** The Rmd above has not been modernized for this
pipeline, and nothing between it and this file recalculates the CVs, so step 05 stages
the table rather than computing it. See `build.sh`'s TODO.


To refresh:

```bash
cp ../../../../MotrpacHumanPreSuspensionAnalysis/data/METABOLOMICS_CVS.rda .
```
