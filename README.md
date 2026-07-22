# MASLD_sex_immune

Scripts and processed data for a meta-analysis of **sex differences in the hepatic immune landscape of MASLD/MASH**, using public bulk-transcriptomic data.

The pipeline recovers each patient's sex from expression, scores immune cell-type and functional-state signatures per sample, and tests — across five cohorts jointly — which immune readouts are systematically sex-biased once fibrosis stage is accounted for. Each headline finding is required to survive two fibrosis codings, FDR correction, leave-one-cohort-out, an independent deconvolution method, and a disease-free (GTEx) control.

The pipeline is written in R (RStudio / HPC).

---

## Layout

```
MASLD_sex_immune/
├── README.md
├── data/                         # deposited processed data (reproduces the analysis)
│   ├── analysis_matrix.csv       #   per-sample signature scores + covariates
│   ├── signature_defs.json       #   the marker/signature definitions used
│   ├── deconv_xcell.csv          #   xCell deconvolution scores
│   └── deconv_mcp.csv            #   MCP-counter deconvolution scores
└── scripts/
    ├── install_packages.R        # one-time: install all dependencies
    ├── 00_config.R               # shared paths, signatures, helpers  <-- edit paths here
    ├── geo_loaders.R             # shared GEO download/parse helpers
    ├── 01_sex_assignment.R       # expression-based sex calls (XIST vs Y panel)
    ├── 02_build_matrix.R         # singscore signature matrix + covariates
    ├── 03_maineffect.R           # sex main-effect mixed models (the core test)
    ├── 04_gtex_control.R         # disease-free GTEx liver control
    ├── 05_export_expr.R          # export linear-scale matrices for deconvolution
    ├── 06_prep_deconv.R          # symbol-level, linear matrices for immunedeconv
    ├── 07_run_deconv.R           # xCell + MCP-counter deconvolution
    ├── 08_concordance.R          # deconvolution vs singscore sex directions
    └── 09_deconfound.R           # BMI/T2D/age deconfounding (within-cohort)
```

## Data

All primary data are public. The four discovery cohorts are on NCBI GEO
(**GSE130970, GSE162694, GSE89632, GSE135251**) and healthy-liver data are from
the **GTEx v8** release. Raw expression is not redistributed here — the scripts
fetch it from those sources. `data/` holds the small processed output needed to
reproduce the results without re-downloading: the per-sample signature matrix,
the signature definitions, and the deconvolution score tables. No controlled-access
or patient-identifiable data are used.

## How to run

Install dependencies once:

```r
source("scripts/install_packages.R")
```

### A. Reproduce the results from the deposited matrix (no downloads)

```r
setwd("scripts")
dir.create("../results", showWarnings = FALSE)
file.copy(list.files("../data", full.names = TRUE), "../results", overwrite = TRUE)
Sys.setenv(MASLD_REALDIR = normalizePath("../results"))
source("03_maineffect.R")
source("08_concordance.R")
source("09_deconfound.R")
```

### B. Full pipeline from raw data

Regenerates everything from GEO/GTEx (downloads several GB the first time).
Edit the paths in `00_config.R` first, then run `01` through `09` in order:

```r
setwd("scripts")
source("01_sex_assignment.R"); source("02_build_matrix.R"); source("03_maineffect.R")
source("04_gtex_control.R");   source("05_export_expr.R");  source("06_prep_deconv.R")
source("07_run_deconv.R");     source("08_concordance.R");  source("09_deconfound.R")
```

## Code availability

Analysis code and the deposited processed data are archived in this repository
(a versioned release is assigned a Zenodo DOI on acceptance).
