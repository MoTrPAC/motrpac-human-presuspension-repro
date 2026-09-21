# GMT builder raw sources

Raw inputs for the buildable GMT leaf data objects. Each is consumed by the
matching `../GMT_*.R` builder to produce a `*.gmt.gz` in `../../data/`. Vendored
copies of the sources in `MotrpacHumanPreSuspensionAnalysis/data-raw/gmt_processing/`.
Provenance and dates below are taken from the original builder-script comments.

| Source file | Provenance | Downloaded / version | Builder | License |
|-------------|-----------|----------------------|---------|---------|
| `CellMarker_2024.txt.gz` | Enrichr — https://maayanlab.cloud/Enrichr/#libraries (CellMarker 2024) | downloaded **2025-01-23** | `GMT_CellMarker.R` | free (cite CellMarker / Enrichr) |
| `Human.MitoCarta3.0.xls` | Broad MitoCarta 3.0 — https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways | MitoCarta **3.0** (no download date recorded) | `GMT_MitoCarta.R` | free (cite Broad MitoCarta3.0) |
| `Kinase_Substrate_Dataset.gz` | PhosphoSitePlus — https://www.phosphosite.org/staticDownloads.action | **v6.7.1.1**; source "Last Modified Fri Nov 17 08:50:20 EST 2023" | `GMT_PSP_kinase.R` | **non-commercial use only**; registration required |
| `ptm.sig.db.all.flanking.human.v2.0.0.gmt.gz` | PTMsigDB — https://proteomics.broadapps.org/ptmsigdb/ | **v2.0.0** (no download date recorded) | `GMT_PTMSigDB.R` | free (cite PTMsigDB v2.0) |

Notes:
- Dates that read "no download date recorded" had no timestamp in the upstream
  builder comments; only the version/source is known.
- To refresh a source, re-download from the URL above and replace the file in
  place; the builder is deterministic given the source.
