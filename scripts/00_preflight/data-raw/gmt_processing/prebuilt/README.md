# Frozen / prebuilt GMT data objects

These `*.gmt.gz` files are **used directly** — they are not rebuilt locally,
because their sources are either login-gated downloads (MSigDB) or a live
web-form query that was snapshotted at a fixed date (RefMet). They are the
data objects themselves, vendored so the pipeline is self-contained.
Provenance/dates below come from the builder-script comments and `scripts/00_preflight/config/gmt_files.tsv`.

| File | Provenance | Version / accessed | Rebuildable here? | License |
|------|-----------|--------------------|-------------------|---------|
| `c2.all.v2023.2.Hs.symbols.gmt.gz` | MSigDB C2 — https://www.gsea-msigdb.org/gsea/msigdb | **v2023.2** | No — download-only (free account/login required) | MSigDB terms |
| `c5.go.v2023.2.Hs.symbols.gmt.gz` | MSigDB C5/GO — https://www.gsea-msigdb.org/gsea/msigdb | **v2023.2** | No — download-only (free account/login required) | MSigDB terms |
| `metabolomics.workbench.refmet.2024.08.07.metabolites.gmt.gz` | RefMet subclass form — https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmetF_form.php | live query **accessed 2024-08-07** | No — needs `HUMAN_FEATURE_TO_GENE` (upstream, non-leaf) + a manual RefMet form submission; see `../GMT_RefMet_subclass.R` | free (stamped snapshot) |

Note: RefMet is listed under the `gmt_files` node but is **not** truly leaf —
its builder depends on the `HUMAN_FEATURE_TO_GENE` object and a manual API
snapshot, so the frozen 2024-08-07 GMT is vendored here instead of rebuilt.
