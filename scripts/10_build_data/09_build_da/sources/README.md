# Step 09 (build_da) raw sources

Raw inputs for the Stage 1 DA stems in `../` that are **not** regenerated inside this
pipeline — external artifacts vendored here so the stems that consume them run offline.
Same convention as `../../06_generate_qc_norm/sources/`.

Everything else this step vendors is *code* (the engine inside `da_common.R`,
`filter_paired_n.R`) and sits alongside the stems rather than here.

## `methylcap_da/` — for `../generate_methylcap_da.R`

MethylCap-seq differential-analysis tables, one file per tissue:
`human-precovid-sed-adu_{t03-edta,t06-muscle,t11-adipose}_epigen-methylcap-seq_da_malax-glmm-acute_v1.2.txt`
(~0.8–1.2 GB each, ~2.7 GB total).

**Provenance:** produced by **Yongchao Ge (yongchao.ge@mssm.edu)**. MethylCap-seq
differential analysis is fit with a **MALAX GLMM** in an external, assay-specific
pipeline; **that modeling code is not included in this repo**, and the dream/voom engine
every other ome shares cannot reproduce it. These tables are therefore the pipeline's
*input*, not its output: `generate_methylcap_da.R` copies them into
`staging/freeze/epigenomics/da/` verbatim and does not modify or re-parse them.

They keep their upstream `malax-glmm-acute_v1.2` filename rather than being renamed to
this pipeline's `dream-acute` token — the model is not dream, and the qc-norm tier vendors
its `v1.2` beta-values under the same reasoning.

**Schema differs from the dream-fit omes.** Columns shared with `.convert_dream_output()`:
`assay`, `feature_id`, `t`, `AveExpr`, `p_value`, `adj_p_value`, `contrast`, `full_model`.
The effect size is **`methylation_diff`, not `logFC`**, and dream's `z.std`,
`degrees_of_freedom` and `logLik` are absent. Nothing is renamed on the way in.

Step 10 (`DA_ASSEMBLE`) is now built, and it **excludes epigen** — matching upstream's
`epigen = FALSE`, since the epigen DA tier is fetched from AWS at call time rather than
shipped as package `.rda`. So these tables still reach no consumer inside this repo, and the
schema clash above is not currently exercised. It is respected structurally all the same:
`../../10_build_da_assemble/differential_analysis_results.R` orders and keys columns through
`intersect()` and never names `logFC`, so a MALAX table would reshape rather than error if
the exclusion filter ever moved.

**Downloaded 2026-07-29** from the production data hub:
`gs://motrpac-data-hub/analysis/human-precovid-sed-adu/c1.3/epigenomics/da/`.
These carry consortium-gated sample/vial labels — do not redistribute outside an
approved data-access agreement.

These files exceed GitHub's 100 MB per-file limit and are gitignored. To fetch or
refresh them:

```bash
make sources GROUPS=methylcap_da
```

`config/copy_from_source.sh` globs on `*methylcap*_da_*.txt`: that bucket folder also holds
the ATAC DA tables, which belong to the separate `atac_da` group below.

## `atac_da/` — for `../generate_atac_da.R`, **only when `RERUN_ATAC=FALSE`**

ATAC-seq differential-analysis tables, one file per tissue:
`human-precovid-sed-adu_{t05-pbmc,t06-muscle}_epigen-atac-seq_da_dream-acute_v1.2.txt`.

Unlike the methylcap tables, these are **not** an external artifact — this pipeline can and
by default does produce them, by fitting dream in `../generate_atac_da.R`. They are
vendored here only to support `RERUN_ATAC=FALSE` in `config/pipeline.env`, which stages the
released tables instead of refitting so that a run can be scoped to the other omes. On that
path `generate_atac_da.R` copies them into the freeze verbatim, does not source
`da_common.R`, and therefore does not depend on steps 05–08.

**Schema is the standard dream output** — same columns `.convert_dream_output()` writes for
every other dream-fit ome, `logFC` included. The methylcap caveat above does not apply.

**Downloaded from** the production data hub:
`gs://motrpac-data-hub/analysis/human-precovid-sed-adu/c1.3/epigenomics/da/`.
These carry consortium-gated sample/vial labels — do not redistribute outside an
approved data-access agreement.

A bare `make sources` skips this group unless `RERUN_ATAC=FALSE`. To fetch it explicitly:

```bash
make sources GROUPS=atac_da
```
