# Processed data

These files are the small processed inputs that ship in the repo so the analysis,
tables and figures reproduce without re-downloading the raw data:

- `analysis_matrix.csv`             – per-sample singscore signatures + covariates
- `signature_defs.json`             – the marker / signature definitions used
- `deconv_xcell.csv`                – xCell deconvolution scores
- `deconv_mcp.csv`                  – MCP-counter deconvolution scores
- `own_cohort_percell_fractions.csv`– in-house CITE-seq per-patient immune fractions (Fig 7)
- `own_cohort_sex.csv`              – patient -> sex map (used by Fig 7 provenance route B)

Raw expression data are public and fetched by the scripts from GEO
(GSE130970, GSE162694, GSE89632, GSE135251), GTEx v8, and the single-cell cohorts.
The in-house fractions are aggregate counts only (no per-cell or identifiable data).
To reproduce the bulk results from these files, see "How to run" in the top-level
README.
