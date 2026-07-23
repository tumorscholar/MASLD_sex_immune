# Processed data

Small processed inputs that ship in the repo. The in-house fractions are
aggregate per-patient counts only (no per-cell or identifiable data):

- `signature_defs.json`              – the marker / signature definitions used
- `own_cohort_percell_fractions.csv` – in-house CITE-seq per-patient immune fractions (Fig 7)
- `own_cohort_sex.csv`               – patient -> sex map (used by Fig 7 provenance route B)

The larger processed tables are **not** committed here; they are regenerated into
the analysis output directory (`MASLD_REALDIR`) when the pipeline runs, and can be
copied back into `data/` before archiving a release if reproduction without
re-downloading raw data is wanted:

- `analysis_matrix.csv` – per-sample singscore signatures + covariates (from `02_build_matrix.R`)
- `deconv_xcell.csv`, `deconv_mcp.csv` – deconvolution scores (from `07_run_deconv.R`)
- the per-donor single-cell fraction CSVs (from the `sc_*` validations)

Raw expression data are public and fetched by the scripts from GEO
(GSE130970, GSE162694, GSE135251, GSE89632, GSE167523, GSE48452), GTEx v8, and the
single-cell cohorts. To reproduce the bulk results, see "Quick start" in the
top-level README.
