# Checks to run against the pipeline output before submission

Open items where a number in the manuscript must be reconciled with what the
analysis actually produces. Verify each when the figures/tables are regenerated
on the HPC, and update the manuscript text to match the pipeline, not the other
way round.

## 1. In-house CITE-seq donor counts (PRIORITY)

The manuscript quotes two different male/female donor counts for the in-house
cohort, and they do not reconcile:

- **Cohort total (Figure 6, Results 3.7):** 19 MASLD donors = **5 male / 14 female**.
- **MAIT functional comparison (Discussion para 2, Supplementary Figure S3):**
  "**seven male and four female** donors", and S3 attributes this to the in-house
  cohort while the Discussion calls it "disease-free liver".

Seven males cannot come from a 5-male cohort. When `sc_build_owncohort.R` /
`sc_mait_functional.R` are re-run, confirm from the actual per-donor table:
  - the true male/female donor split of the in-house cohort (drives Figure 6 and
    Results 3.7), and
  - which donors and which tissue (MASLD in-house vs a disease-free set) go into
    the MAIT per-cell functional comparison (drives Discussion para 2 and S3).

Then make the two statements consistent (either both = the real in-house split,
or the functional comparison is explicitly a different/disease-free set with its
own n). Files touching these numbers: `fig6_singlecell.R`, `sc_mait_functional.R`,
`sc_build_owncohort.R`; manuscript spots: Results 3.7, Discussion para 2,
Fig 6 legend, Supp Fig S3 legend.

## 2. Other live numbers to confirm against the regenerated CSVs

These should already agree, but re-check after the run since they are quoted to
2–3 sig figs in the text:
  - Sex-call concordance 95–98% (per-cohort; Results 3.1, Fig 1B) from
    `per_sample_sex_calls.csv`.
  - Six-cohort meta β / 95% CI / I² for every headline readout (Results 3.2,
    Table 2, Supp Table S2) from `13_meta_random.R` output.
  - GTEx balanced-resampling: receptor MAIT nominal p = 0.049, direction held in
    100/100 subsamples (Results 3.5, Supp Table S3).
  - Single-cell pooled p = 0.10 (mixed model) / 0.12 (van Elteren), Andrews
    p = 0.019, HLiCA p = 0.13, median 50 cells/donor (Results 3.7, Supp Tables S5).
  - Deconvolution concordance 13/18 (72%) (Results 3.6, Table S4).

## 1b. UPDATE from the first full HPC run (2026-07-31)

`fig6_singlecell.R` now prints `[CHECK] in-house cohort donor split: 5 male /
14 female (n=19)` — so the **5 M / 14 F** figure in Results 3.7 and the Fig 6
legend is correct, and it is the **7 M / 4 F** in Discussion para 2 + Supp Fig S3
that must be fixed (or explicitly attributed to a different/disease-free set with
its own n). Reconcile toward 5 M / 14 F unless the functional comparison genuinely
used a different sample set.

## 4. GSE167523 dropped in the run — MUST FIX (blocks the six-cohort claim)

The run logged `!! FAILED GSE167523 : no tabular suppl for GSE167523`. The
per-cohort `try()` lets the step finish "OK", so the SUMMARY hides it, but
GSE167523's count table never downloaded and it was excluded. That silently
drops the study to five cohorts and breaks numbers the manuscript quotes
(n = 616, the "94.9% in GSE167523" concordance in Results 3.1, the six-cohort
meta and I² values). Do not trust any regenerated table/figure until GSE167523
is back in. Diagnose:

```r
source("00_config.R")
list.files(file.path(GEO_CACHE, "GSE167523"), recursive = TRUE)  # what actually downloaded
read.csv(file.path(REALDIR, "analysis_matrix.csv")) |> (\(d) table(d$cohort))()  # is GSE167523 present?
```

If the cache dir is empty, the supplementary download failed (re-fetch, or
download the count table on a login node into `GEO_CACHE/GSE167523/`). If it has
a file the loader didn't recognise (not `.csv/.tsv/.txt/.tab`, e.g. `.xlsx` or an
oddly named counts file), place a tabular copy there so `rnaseq_matrix()` picks
it up, then re-run `01`, `02`, and everything downstream.

## 5. Deconfounding (step 09) failed — GSE89632 metabolic covariates empty

`09_deconfound.R` needs per-patient BMI / T2D / age for GSE89632. It failed;
`09` now stops with a clear message printing the non-NA counts. Check:

```r
d <- read.csv(file.path(REALDIR, "analysis_matrix.csv"))
s <- d[d$cohort == "GSE89632", ]; colSums(!is.na(s[, c("bmi","t2d","age")]))
```

If these are 0, the covariates were not parsed from GEO. Supply a per-patient
metabolic table (merge into `analysis_matrix.csv`) to regenerate the check.
NOTE: Supp Fig S4 showed "OK" only because a `deconf_within_cohort.csv` from a
PRIOR run was still in `MASLD_REALDIR` — the S4 figure is therefore STALE until
`09` runs successfully. Delete the old file if you want S4 to fail loudly rather
than plot old numbers.

## 3. Software versions (author-side, not a figure run)

Fill the 16 `[version]` placeholders in the manuscript's Supplementary methods
"Software" section from `sessionInfo()` captured on the HPC after
`install_packages.R`.
