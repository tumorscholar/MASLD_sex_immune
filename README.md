# MASLD_sex_immune

Scripts and processed data for a standalone meta-analysis of **sex differences in the hepatic immune landscape of MASLD/MASH**, using only public bulk-transcriptomic data. Public data in, published signatures out — the analysis is independent of any single-cell data.

The pipeline recovers each patient's sex from expression, scores immune cell-type and functional-state signatures per sample, and tests — across five cohorts jointly — which immune readouts are systematically sex-biased once fibrosis stage is accounted for. Every headline finding is required to survive two fibrosis codings, FDR correction, leave-one-cohort-out, an independent deconvolution method, and a disease-free (GTEx) control.

The pipeline is written in **R** (RStudio / HPC). For a plain-language walk-through of every method and statistic, see **[`METHODS_EXPLAINED.md`](METHODS_EXPLAINED.md)**.

---

## Layout

```
MASLD_sex_immune/
├── README.md
├── METHODS_EXPLAINED.md          # plain-language guide to the stats
├── data/                         # deposited processed data (small; reproduces the paper)
│   ├── analysis_matrix.csv       #   per-sample signature scores + covariates
│   ├── signature_defs.json       #   the marker/signature definitions used
│   ├── deconv_xcell.csv          #   xCell deconvolution scores
│   └── deconv_mcp.csv            #   MCP-counter deconvolution scores
└── scripts/                      # the R pipeline
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
    ├── 09_deconfound.R           # BMI/T2D/age deconfounding (within-cohort)
    └── 10_verify_vs_python.R     # optional: reproducibility cross-check
```

## Data

All primary data are public. The four discovery cohorts are on NCBI GEO
(**GSE130970, GSE162694, GSE89632, GSE135251**) and healthy-liver data are from
the **GTEx v8** release. Raw expression is not redistributed here — the scripts
fetch it from those sources. What is deposited in `data/` is the small *processed*
output needed to reproduce the paper's statistics without re-downloading anything:
the per-sample signature matrix (`analysis_matrix.csv`), the signature definitions,
and the two deconvolution score tables. No controlled-access or patient-identifiable
data are used.

## How to run

Install dependencies once:

```r
source("scripts/install_packages.R")
```

### A. Reproduce the paper's statistics (no downloads)

Runs the models on the deposited `analysis_matrix.csv` — this reproduces the
manuscript's sex effects, deconvolution concordance and deconfounding results.

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

Regenerates everything from GEO/GTEx (slow — downloads several GB the first time).
Edit the paths in `00_config.R` first, then:

```r
setwd("scripts")
for (s in sprintf("%02d_*.R", 1:9)) NULL   # run 01 through 09 in order:
source("01_sex_assignment.R"); source("02_build_matrix.R"); source("03_maineffect.R")
source("04_gtex_control.R");   source("05_export_expr.R");  source("06_prep_deconv.R")
source("07_run_deconv.R");     source("08_concordance.R");  source("09_deconfound.R")
```

## Reproducibility note

The deposited `analysis_matrix.csv` is the reference input, and running mode A on
it reproduces the manuscript numbers. The full fetch-and-build pipeline (mode B)
regenerates that matrix from the raw data; because gene-ID mapping libraries differ
in coverage, the rebuilt singscore magnitudes for the RNA-seq cohorts can differ
marginally from the deposited matrix (worst case ≈0.1 SD on a pooled sex
coefficient). This does not change any effect direction or conclusion — every
sex direction and the GTEx, deconvolution and deconfounding results reproduce
exactly. `10_verify_vs_python.R` documents this cross-check. The analysis was
originally developed in Python; that implementation is retained by the authors and
is available on request.

## Citation / code availability

Analysis code and the deposited processed data are archived in this repository
(a versioned release is assigned a Zenodo DOI on acceptance). The in-house CITE-seq
validation cohort is deposited separately at NCBI GEO (expression) and the European
Genome-phenome Archive (controlled-access sequence data); accession numbers are in
the manuscript.
