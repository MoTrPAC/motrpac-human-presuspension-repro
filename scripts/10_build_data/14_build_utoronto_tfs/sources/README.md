# Step 14 (build_utoronto_tfs) raw sources

Same convention as `../../05_stage_metabolomics_cvs/sources/` and
`../../09_build_da/sources/`: inputs this pipeline does not regenerate, vendored so the
step runs without reaching outside the repo.

## `DatabaseExtract_v_1.01.txt`

The Human Transcription Factors database extract — one row per candidate human TF, with
the curated `Is TF?` call and its supporting evidence columns.

| | |
|---|---|
| **Source** | <https://humantfs.ccbr.utoronto.ca/download.php> ("Current lists" → the full database extract) |
| **Citation** | Lambert SA, Jolma A, Campitelli LF, et al. *The Human Transcription Factors.* Cell. 2018;172(4):650-665. |
| **Version** | v1.01 (the version is in the filename; the site publishes no separate release date) |
| **Vendored** | 2026-08-11 |
| **md5** | `43c81af9799cc41e5a95bb6d810fe9b0` |
| **Shape** | 2,765 data rows x 29 columns — column 1 is an unnamed row index, dropped on read, which is why the shipped object is 2,765 x 28 |
| **License** | Academic use, per the download page. Not consortium data — this is a public external reference, unlike the gated sources under `../../09_build_da/sources/`. |

The first column is an **unnamed row index**, which is why `../UTORONTO_TFs.R` drops it
(`[-1]`) — the same thing the upstream generator does. `check.names = TRUE` is likewise
deliberate and load-bearing: it is what turns `Is TF?` into the `Is.TF.` the filter names
and `HGNC symbol` into `HGNC.symbol`.

To refresh, download the current extract from the link above, drop it here under its
versioned filename, and update the version/date/md5 rows in this table.
