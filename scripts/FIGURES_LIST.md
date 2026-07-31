# Figure and table catalogue — sex × MASLD hepatic-immune study

Every figure and table is regenerated from the pipeline output CSVs. Outputs are
written under your analysis output directory (`MASLD_REALDIR`):

- **figures** → `MASLD_REALDIR/figures/` as **PDF (vector) + TIFF (600 dpi, LZW)**
  (JHEP/Elsevier format: RGB, Helvetica/Arial, no PNG).
- **tables** → `MASLD_REALDIR/tables/` as CSV plus a combined `MASLD_sex_tables.xlsx`.

One self-contained figure script per Results section: open it in RStudio and
Source it (each reads its own CSVs and writes its own PDF/TIFF; no need to source
`run_all.R` first). `run_all.R` runs them all in order.

## Main figures (6) — one per Results section

| Fig | Section | Script | Panels | Source data |
|-----|---------|--------|--------|-------------|
| **1** | 3.1 | `fig1_sex_assignment.R`   | A calls, B validation, C bimodal | `per_sample_sex_calls.csv` |
| **2** | 3.2 | `fig2_main_effects.R`     | A forest (all readouts), B distributions (headline) | `maineffect_results.csv`, `maineffect_concordance.csv`, `analysis_matrix.csv` |
| **3** | 3.3 | `fig3_mait_identity.R`    | A receptor specificity, B effect size by marker set | `analysis_matrix.csv`, `maineffect_results.csv` |
| **4** | 3.5 | `fig4_gtex_specificity.R` | A MASLD-vs-GTEx, B cross-cohort robustness | `gtex_healthy_sex.csv`, `maineffect_results.csv` |
| **5** | 3.6 | `fig5_deconvolution.R`    | single panel | `deconv_concordance.csv` |
| **6** | 3.7 | `fig6_singlecell.R`       | one panel per readout | in-house `*_per_donor_fractions.csv` (`OWN_CITESEQ_CSV`) |

## Main tables (2)

| Table | Script | Source data |
|-------|--------|-------------|
| **1** Datasets / cohort characteristics | `make_tables.R` | `analysis_matrix.csv` |
| **2** Headline sex main effects | `make_tables.R` | `maineffect_results.csv`, `maineffect_concordance.csv` |

## Supplementary figures (S1–S8)

| SFig | Section | Script |
|------|---------|--------|
| **S1** Analysis workflow schematic | Methods 2.1 | `make_workflow_schematic.R` |
| **S2** Bulk random-effects meta forest | 3.2 | `13_meta_random.R` |
| **S3** MAIT functional state by sex | 3.3 | `sc_mait_functional.R` |
| **S4** Metabolic deconfounding | 3.4 | `suppfig_deconfounding.R` |
| **S5** Per-cohort single-cell analysis | 3.7 | `11_supp_figures.R` |
| **S6** Single-cell per-cohort effect sizes | 3.7 | `11_supp_figures.R` (from `sc_meta_forest.R` / `16_sc_pooled_mait.R`) |
| **S7** Guilliams capture QC + sex×diet check | 3.7 | `11_supp_figures.R` |
| **S8** Single-cell random-effects meta forest | 3.7 | `13_meta_random.R` |

## Supplementary tables (S1–S5)

| STable | Section | Script |
|--------|---------|--------|
| **S1** Marker genes for each signature | Methods 2.3 | `make_tables.R` |
| **S2** Random-effects meta β, 95% CI, I² | 3.2 | `13_meta_random.R` → `make_tables.R` |
| **S3** Disease-specificity (MASLD vs GTEx) | 3.5 | `make_tables.R` |
| **S4** Benchmarked deconvolution concordance | 3.6 | `make_tables.R` |
| **S5** Single-cell validation cohorts | 3.7 | `make_tables.R` |

## Notes
- §3.4 (metabolic deconfounding) has no main figure; it is Supplementary Figure S4.
- Figure 4 panels: A = MASLD-vs-GTEx disease-specificity, B = cross-cohort robustness.
- Supplementary figures are numbered by order of first citation (workflow S1 first).
- The supplementary figure scripts (`11_supp_figures.R`, `sc_mait_functional.R`)
  still emit their original filenames; only the new main scripts and
  `suppfig_deconfounding.R` (S4) / `make_workflow_schematic.R` (S1) self-label to
  the final scheme.
