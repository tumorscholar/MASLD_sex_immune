# Figure and table catalogue – sex x MASLD hepatic-immune study

Every figure and table in the manuscript is regenerated from the pipeline output
CSVs. Outputs are written under your analysis output directory (`MASLD_REALDIR`):
figures to the `figures/` subfolder as **PDF (vector) + TIFF (300 dpi, LZW) + PNG**,
and tables to the `tables/` subfolder as CSV plus a combined `MASLD_sex_tables.xlsx`.

## How to generate everything

```r
setwd("scripts")
source("run_all.R")     # analysis -> tables -> figures, in order
```

or piece by piece:

```r
source("make_tables.R")        # Tables 1-5 (+ 5b)
source("10_figures.R")         # Fig 1-6 + the in-house panel of Supp Fig S1
source("sc_meta_dotplot.R")    # Fig 7 (combined in-house + public single-cell)
source("11_supp_figures.R")    # Supp Fig S1 (public panels), S2, S3
source("sc_mait_functional.R") # Supp Fig S4 (MAIT functional state), needs the Andrews object
```

Set `MASLD_REALDIR` if your paths differ from `00_config.R`. Fig 7 and the in-house
S1 panel find `data/own_cohort_percell_fractions.csv` automatically.

## Main figures

| ID | Title | Script | Source data |
|----|-------|--------|-------------|
| **Fig 1** | Expression-based sex assignment | `10_figures.R` | `per_sample_sex_calls.csv` |
| **Fig 2** | Pooled sex main-effect forest | `10_figures.R` | `maineffect_results.csv`, `maineffect_concordance.csv` |
| **Fig 3** | Per-sample distributions of headline readouts | `10_figures.R` | `analysis_matrix.csv`, `maineffect_results.csv` |
| **Fig 4** | MAIT receptor-identity specificity | `10_figures.R` | `analysis_matrix.csv`, `maineffect_results.csv` |
| **Fig 5** | Disease-specificity + robustness (MASLD vs GTEx) | `10_figures.R` | `gtex_healthy_sex.csv`, `deconf_within_cohort.csv`, `maineffect_results.csv` |
| **Fig 6** | Deconvolution concordance | `10_figures.R` | `deconv_concordance.csv` |
| **Fig 7** | Single-cell validation, in-house + four public cohorts (one combined dot plot) | `sc_meta_dotplot.R` | in-house + per-cohort `*_per_donor_fractions.csv` |

## Supplementary figures

| ID | Title | Script | Source data |
|----|-------|--------|-------------|
| **S1** | Per-cohort single-cell validation (in-house + 4 public, one panel each) | in-house panel: `10_figures.R`; public panels: `11_supp_figures.R` | per-cohort fractions |
| **S2** | Single-cell meta-analysis forest (Cliff's δ per cohort + pooled) | `11_supp_figures.R` (from `sc_meta_forest.R` output) | `meta/sc_meta_effects.csv` |
| **S3** | Guilliams capture QC + sex x diet confound check | `11_supp_figures.R` | `guilliams/out/Guilliams_clean_per_patient.csv` |
| **S4** | MAIT functional-state scores by sex (abundance not phenotype) | `sc_mait_functional.R` | Andrews healthy MAIT pseudobulk |

## Tables

| ID | Title | Script | Source data |
|----|-------|--------|-------------|
| **Table 1** | Datasets / cohort characteristics | `make_tables.R` | `analysis_matrix.csv` |
| **Table 2** | Headline sex main effects | `make_tables.R` | `maineffect_results.csv`, `maineffect_concordance.csv` |
| **Table 3** | Disease-specificity (MASLD vs GTEx) | `make_tables.R` | `gtex_healthy_sex.csv` |
| **Table 4** | Benchmarked deconvolution concordance | `make_tables.R` | `deconv_concordance.csv` |
| **Table 5** | Single-cell validation summary (cohort x readout) | `make_tables.R` | per-cohort fractions + `own_cohort_percell_fractions.csv` |
| **Table 5b** | Single-cell cohort donor counts (N, M/F) for §3.7 / Fig 7 | `make_tables.R` | derived from the per-cohort fractions |

## Notes

- Fig 7 is the single combined single-cell figure (in-house cohort as the first
  column, then the four public cohorts). `10_figures.R` also emits the in-house
  cohort on its own as the in-house panel of Supp Fig S1.
- Every figure is produced identically in PDF / TIFF / PNG, so any row can move
  between main and supplementary freely.
