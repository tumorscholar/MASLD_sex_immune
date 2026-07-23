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
  "remotes",      # install tools from GitHub
  "Matrix",       # sparse matrices (single-cell)
  "Seurat",       # single-cell objects and processing
  "msigdbr",      # gene-set collections for fgsea
  "openxlsx",     # combined tables workbook (make_tables.R)
  "metafor"       # random-effects meta-analysis (13_meta_random.R)
)
for (p in cran)
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

## Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc <- c(
  "GEOquery",     # download + parse GEO series
  "org.Hs.eg.db", # Entrez/Ensembl <-> HGNC symbol mapping
  "AnnotationDbi",
  "edgeR",        # pseudobulk differential expression (sc_mait_functional.R)
  "limma",        # edgeR dependency
  "fgsea",        # gene-set enrichment (sc_mait_functional.R)
  "SingleCellExperiment", # reference container for MuSiC (12_liver_deconv.R)
  "TOAST"         # MuSiC dependency (must be installed before MuSiC)
)
for (p in bioc)
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)

## Tools installed from GitHub
if (!requireNamespace("xCell", quietly = TRUE))          # deconvolution (07_run_deconv.R)
  remotes::install_github("dviraran/xCell", upgrade = "never")
if (!requireNamespace("MCPcounter", quietly = TRUE))     # deconvolution (07_run_deconv.R)
  try(remotes::install_github("ebecht/MCPcounter", subdir = "Source", upgrade = "never"))
if (!requireNamespace("schard", quietly = TRUE))         # read h5ad (sc_validate_hlica.R)
  try(remotes::install_github("cellgeni/schard", upgrade = "never"))
if (!requireNamespace("MuSiC", quietly = TRUE))          # liver-reference deconvolution (12_liver_deconv.R)
  try(remotes::install_github("xuranw/MuSiC", upgrade = "never"))

message("install_packages.R done.")
