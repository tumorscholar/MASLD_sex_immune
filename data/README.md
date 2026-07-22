# Processed data

These files are the processed output of the pipeline and let anyone reproduce the
analysis without re-downloading the raw data:

- `analysis_matrix.csv`  — per-sample singscore signatures + covariates
- `signature_defs.json`  — the marker / signature definitions used
- `deconv_xcell.csv`     — xCell deconvolution scores
- `deconv_mcp.csv`       — MCP-counter deconvolution scores

Raw expression data are public and fetched by the scripts from GEO
(GSE130970, GSE162694, GSE89632, GSE135251) and GTEx v8. To reproduce the results
from these files, see "How to run → A" in the top-level README.
