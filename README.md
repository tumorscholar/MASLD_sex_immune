# MASLD_sex_immune

Scripts and processed data for a meta-analysis of **sex differences in the hepatic immune landscape of MASLD/MASH**.

The study is built entirely on public data in, published signatures out. The bulk pipeline recovers each patient's sex from expression, scores immune cell-type and functional-state signatures per sample, and tests across four MASLD cohorts jointly which immune readouts are sex-biased once fibrosis stage is accounted for. Each headline finding is then required to survive two fibrosis codings, FDR correction, leave-one-cohort-out, an independent deconvolution method, a metabolic (BMI/T2D/age) deconfounding step, and a disease-free (GTEx) control. Finally the directions are validated at single-cell resolution in our own in-house CITE-seq cohort and four independent public cohorts (section 3.7; combined in Fig 7).

Everything is written in R and runs in RStudio on the HPC or with plain `Rscript`.

---

## Two axes of the result

1. A **baseline male MAIT axis** – higher MAIT abundance in males, present already in disease-free liver (GTEx), receptor-identity specific (SLC4A10/TRAV1-2), and reproduced across bulk, deconvolution, GTEx, and every single-cell cohort tested.
2. A **disease-emergent female programme** – Treg / CD8 / cDC / Th1 signatures higher in females in MASLD but not in healthy liver, i.e. they appear with disease rather than being baseline.

---

## Layout

```
MASLD_sex_immune/
├── README.md
├── run_all.R  ->  scripts/run_all.R          # end-to-end entry point
├── data/                                      # small processed data (ships in the repo)
│   ├── analysis_matrix.csv                    #   per-sample signature scores + covariates
│   ├── signature_defs.json                    #   marker / signature definitions
│   ├── deconv_xcell.csv, deconv_mcp.csv       #   deconvolution scores
│   ├── own_cohort_percell_fractions.csv       #   in-house CITE-seq per-patient fractions (Fig 7)
│   └── own_cohort_sex.csv                      #   patient -> sex map (Fig 7 provenance route B)
├── results_R/                                 # regenerated outputs (git-ignored except this note)
│   ├── tables/                                #   Table 1-4 CSVs + xlsx
│   └── figures/                               #   Fig 1-8 + Supp as PDF / TIFF / PNG
└── scripts/
    ├── run_all.R                 # runs the whole pipeline in order (analysis -> tables -> figures)
    ├── install_packages.R        # one-time: install all dependencies
    ├── 00_config.R               # shared paths, signatures, helpers – edit paths here
    ├── geo_loaders.R             # shared GEO download / parse helpers
    ├── 01_sex_assignment.R       # expression-based sex calls (XIST vs Y panel)
    ├── 02_build_matrix.R         # singscore signature matrix + covariates
    ├── 03_maineffect.R           # sex main-effect mixed models (the core test)
    ├── 04_gtex_control.R         # disease-free GTEx liver control (3 models: age / covariate / balanced)
    ├── 05_export_expr.R          # export linear-scale matrices for deconvolution
    ├── 06_prep_deconv.R          # symbol-level linear matrices for immunedeconv
    ├── 07_run_deconv.R           # xCell + MCP-counter deconvolution
    ├── 08_concordance.R          # deconvolution vs singscore sex directions
    ├── 09_deconfound.R           # BMI / T2D / age deconfounding (within cohort)
    ├── sc_build_*.R              # build one Seurat object per public single-cell cohort
    ├── sc_validate_*.R           # gate cell types + test sex directions per cohort
    ├── sc_guilliams_clean.R      # patient-level re-analysis (removes a capture artifact)
    ├── sc_mait_functional.R      # male-vs-female MAIT functional state (abundance not phenotype)
    ├── sc_build_owncohort.R      # in-house CITE-seq fractions (Fig 7 provenance)
    ├── sc_meta_forest.R          # Cliff's delta + inverse-variance pooled single-cell meta
    ├── make_tables.R             # Table 1-4 from the pipeline CSVs
    ├── 10_figures.R              # Fig 1-7 (PDF + TIFF + PNG)
    ├── sc_meta_dotplot.R         # Fig 7 (combined in-house + public single-cell)
    ├── 11_supp_figures.R         # Supp Fig S1-S3
    └── FIGURES_LIST.md           # figure + table catalogue (main vs supplementary)
```

## Data

All primary data are public. The four MASLD cohorts are on NCBI GEO
(**GSE130970, GSE162694, GSE89632, GSE135251**); healthy-liver data are from
**GTEx v8**; the single-cell validation cohorts are Andrews 2024 (GSE243981),
HLiCA, Guilliams 2022 and Ramachandran 2019. Raw data are not redistributed here;
the scripts fetch them from source. `data/` holds only the small processed tables
needed to reproduce the results without re-downloading. No controlled-access or
patient-identifiable data are used. The in-house CITE-seq per-patient fractions
(`data/own_cohort_percell_fractions.csv`) are aggregate counts only.

## How to run

Install dependencies once:

```r
source("scripts/install_packages.R")
```

### Everything, in order

From the repo root:

```r
setwd("scripts")
source("run_all.R")
```

`run_all.R` runs analysis -> tables -> figures and prints a pass/fail/skip summary.
Stages are switched with environment variables (all default sensibly):

| variable | default | what it controls |
|----------|---------|------------------|
| `RUN_BULK`    | on  | bulk pipeline `01`-`09` |
| `RUN_SC`      | off | single-cell validations (need the large raw downloads) |
| `RUN_TABLES`  | on  | `make_tables.R` |
| `RUN_FIGURES` | on  | `10_figures.R`, `sc_meta_dotplot.R`, `11_supp_figures.R` |

So a quick "tables + figures from what is already computed" run is just the default
with `RUN_BULK=0 RUN_SC=0`.

### Reproduce the bulk results from the deposited matrix (no downloads)

```r
setwd("scripts")
dir.create("../results_R", showWarnings = FALSE)
file.copy(list.files("../data", full.names = TRUE), "../results_R", overwrite = TRUE)
Sys.setenv(MASLD_REALDIR = normalizePath("../results_R"))
source("03_maineffect.R"); source("08_concordance.R"); source("09_deconfound.R")
source("make_tables.R")
```

### Full pipeline from raw data

Edit the paths in `00_config.R`, then either run `run_all.R` with `RUN_SC=1`, or
run the numbered scripts `01`-`09`, the `sc_*` scripts, `make_tables.R`, and the
figure scripts in that order.

## Figure 7 (single-cell validation) provenance

Figure 7 is a single combined dot plot: the in-house CITE-seq cohort as the first
column (section 3.7, validation in our own data) followed by the four public
single-cell cohorts (validation in external data). It is produced by
`sc_meta_dotplot.R`. The in-house column is drawn from
`data/own_cohort_percell_fractions.csv`; `sc_build_owncohort.R` regenerates that
table from the full annotated CITE-seq Seurat object (route A, the version used in
the paper; RNA + protein cell-type labels) and offers a self-contained marker-panel
cross-check (route B). The MAIT / Treg / CD8 / cDC gates are identical across the
in-house and public cohorts, so every column of Figure 7 is methodologically matched.

## Code availability

Analysis code and the deposited processed data are archived in this repository.
A versioned release is assigned a Zenodo DOI on acceptance.
