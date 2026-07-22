# Deposited processed data

These small files are the processed output of the pipeline and are what let anyone
reproduce the manuscript's statistics without re-downloading the raw data. Place the
following four files here (copied from your analysis working directory on the HPC):

- `analysis_matrix.csv`  — per-sample singscore signatures + covariates (from 02_build_matrix.R)
- `signature_defs.json`  — the marker / signature definitions used (from 02_build_matrix.R)
- `deconv_xcell.csv`     — xCell deconvolution scores (from 07_run_deconv.R)
- `deconv_mcp.csv`       — MCP-counter deconvolution scores (from 07_run_deconv.R)

Raw expression data are NOT stored here — they are public and fetched by the scripts
from GEO (GSE130970, GSE162694, GSE89632, GSE135251) and GTEx v8.

To reproduce the paper's numbers from these files, see "How to run → A" in the top-level README.
