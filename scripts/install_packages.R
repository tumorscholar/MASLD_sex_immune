## install_packages.R -------------------------------------------------------
## One-time setup. Run this ONCE in your RStudio session on the HPC before
## running the pipeline. Safe to re-run (skips anything already installed).
## ---------------------------------------------------------------------------
options(repos = "https://cloud.r-project.org")

cran <- c(
  "data.table",   # fast table I/O (large GTEx / expression matrices)
  "lme4",         # linear mixed-effects models
  "lmerTest",     # p-values for lmer fixed effects (Satterthwaite)
  "ggplot2",      # figures
  "jsonlite",     # write signature_defs.json
  "remotes"       # install the deconvolution tools from GitHub
)
for (p in cran)
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

## Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc <- c(
  "GEOquery",     # download + parse GEO series (the R equivalent of GEOparse)
  "org.Hs.eg.db", # Entrez/Ensembl <-> HGNC symbol mapping (replaces mygene)
  "AnnotationDbi"
)
for (p in bioc)
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)

## Deconvolution tools (used by 07_run_deconv.R)
if (!requireNamespace("xCell", quietly = TRUE))
  remotes::install_github("dviraran/xCell", upgrade = "never")
if (!requireNamespace("MCPcounter", quietly = TRUE))
  try(remotes::install_github("ebecht/MCPcounter", subdir = "Source", upgrade = "never"))

message("install_packages.R done.")
