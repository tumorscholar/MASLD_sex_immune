# Sex differences in the hepatic immune landscape of MASLD

Analysis code and processed data for a transcriptomic meta-analysis of sex
differences in the hepatic immune compartment of metabolic dysfunction-associated
steatotic liver disease (MASLD / MASH), built entirely on public data.

![Language](https://img.shields.io/badge/language-R-276DC3)
![License](https://img.shields.io/badge/license-MIT-green)

## Overview

The pipeline assigns each sample's sex from expression, scores immune cell-type and
functional-state signatures per sample, and estimates the sex effect for every
readout across four MASLD cohorts jointly, adjusting for fibrosis stage. Each
headline result is required to survive two fibrosis codings, false-discovery-rate
correction, leave-one-cohort-out resampling, an independent deconvolution method,
a metabolic (BMI / type-2 diabetes / age) deconfounding step, and a disease-free
(GTEx) control. The sex directions are then validated at single-cell resolution in
an in-house CITE-seq cohort and four independent public single-cell cohorts.

The analysis resolves two axes of sex difference: a **baseline male MAIT-cell bias**,
present already in disease-free liver and specific to the receptor-defined
(SLC4A10 / TRAV1-2) MAIT programme; and a **disease-emergent, female-skewed
programme** (regulatory T cells, CD8 T cells and conventional dendritic cells) that
appears with disease rather than at baseline.

## Repository structure

```
MASLD_sex_immune/
├── README.md
├── LICENSE
├── data/                      Processed inputs that ship with the repository
│   ├── analysis_matrix.csv        per-sample signature scores and covariates
│   ├── signature_defs.json        marker and signature definitions
│   ├── deconv_xcell.csv           xCell deconvolution scores
│   ├── deconv_mcp.csv             MCP-counter deconvolution scores
│   ├── own_cohort_percell_fractions.csv   in-house CITE-seq per-patient fractions
│   └── own_cohort_sex.csv         in-house patient-to-sex map
├── results_R/                 Regenerated outputs (figures and tables; git-ignored)
└── scripts/
    ├── run_all.R                  end-to-end runner (analysis, tables, figures)
    ├── install_packages.R         install all dependencies
    ├── 00_config.R                shared paths, signatures and helpers (edit paths here)
    ├── geo_loaders.R              GEO download and parsing helpers
    ├── 01_sex_assignment.R        expression-based sex calls (XIST and Y-panel)
    ├── 02_build_matrix.R          singscore signature matrix and covariates
    ├── 03_maineffect.R            sex main-effect mixed models
    ├── 04_gtex_control.R          disease-free GTEx liver control
    ├── 05_export_expr.R           linear-scale matrices for deconvolution
    ├── 06_prep_deconv.R           symbol-level matrices for immunedeconv
    ├── 07_run_deconv.R            xCell and MCP-counter deconvolution
    ├── 08_concordance.R           deconvolution vs signature concordance
    ├── 09_deconfound.R            BMI / type-2 diabetes / age deconfounding
    ├── sc_build_andrews.R         build a Seurat object (Andrews 2024)
    ├── sc_build_guilliams.R       build a Seurat object (Guilliams / Liver Cell Atlas)
    ├── sc_build_ramachandran.R    build a Seurat object (Ramachandran 2019)
    ├── sc_validate_andrews.R      gate cell types and test sex directions (Andrews)
    ├── sc_validate_guilliams.R    gate cell types and test sex directions (Guilliams)
    ├── sc_validate_hlica.R        gate cell types and test sex directions (HLiCA)
    ├── sc_validate_ramachandran.R gate cell types and test sex directions (Ramachandran)
    ├── sc_guilliams_clean.R       patient-level re-analysis of the Guilliams cohort
    ├── sc_build_owncohort.R       in-house CITE-seq per-patient fractions
    ├── sc_mait_functional.R       male vs female MAIT functional-state comparison
    ├── sc_meta_forest.R           Cliff's delta and inverse-variance single-cell meta-analysis
    ├── make_tables.R              Tables 1 to 5
    ├── 10_figures.R               Figures 1 to 6
    ├── sc_meta_dotplot.R          Figure 7 (combined single-cell validation)
    ├── 11_supp_figures.R          Supplementary Figures S1 to S3
    └── FIGURES_LIST.md            figure and table catalogue
```

## Requirements

R (>= 4.1) with the packages listed in `scripts/install_packages.R`, including
`data.table`, `lme4`, `lmerTest`, `singscore`, `immunedeconv` (xCell, MCP-counter),
`Seurat`, `ggplot2`, `fgsea`, `msigdbr` and `openxlsx`. Install them once with:

```r
source("scripts/install_packages.R")
```

## Quick start

Run everything in dependency order (analysis, then tables, then figures) from the
repository root:

```r
setwd("scripts")
source("run_all.R")
```

`run_all.R` prints a pass / fail / skip summary. Stages are controlled by environment
variables, all with sensible defaults:

| Variable      | Default | Controls |
|---------------|---------|----------|
| `RUN_BULK`    | on      | bulk pipeline (`01`–`09`) |
| `RUN_SC`      | off     | single-cell validations (require the large raw downloads) |
| `RUN_TABLES`  | on      | `make_tables.R` |
| `RUN_FIGURES` | on      | `10_figures.R`, `sc_meta_dotplot.R`, `11_supp_figures.R` |

### Reproduce the bulk results without downloading raw data

The deposited `data/` files reproduce the core results directly:

```r
setwd("scripts")
dir.create("../results_R", showWarnings = FALSE)
file.copy(list.files("../data", full.names = TRUE), "../results_R", overwrite = TRUE)
Sys.setenv(MASLD_REALDIR = normalizePath("../results_R"))
source("03_maineffect.R"); source("08_concordance.R"); source("09_deconfound.R")
source("make_tables.R")
```

### Full pipeline from raw data

Edit the paths in `00_config.R`, then run `run_all.R` with `RUN_SC=1`, or run the
numbered scripts `01`–`09`, the `sc_*` scripts, `make_tables.R` and the figure
scripts in that order. `FIGURES_LIST.md` maps every figure and table to the script
that produces it.

## Data availability

All primary data are public. Raw data are not redistributed here; the scripts
retrieve them from source, and `data/` holds only the small processed inputs needed
to reproduce the results. No controlled-access or patient-identifiable data are used;
the in-house CITE-seq file contains aggregate per-patient counts only.

| Source | Accession |
|--------|-----------|
| MASLD bulk cohorts | GEO GSE130970, GSE162694, GSE89632, GSE135251 |
| Disease-free liver | GTEx v8 |
| Single-cell (Andrews 2024) | GEO GSE243981 |
| Single-cell (HLiCA, Guilliams 2022, Ramachandran 2019) | published atlases |

## Citation

If you use this code, please cite the associated manuscript (in preparation). A
versioned release will be archived with a Zenodo DOI on acceptance.

## License

Released under the MIT License. See [LICENSE](LICENSE).
