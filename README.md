# Sex differences in the hepatic immune landscape of MASLD

Analysis code for a transcriptomic meta-analysis of sex differences in the hepatic
immune compartment of metabolic dysfunction-associated steatotic liver disease
(MASLD) and its inflammatory form, metabolic dysfunction-associated steatohepatitis
(MASH) — the current nomenclature for NAFLD/NASH and the most common chronic liver
disease worldwide. MASLD is markedly sex-dimorphic in prevalence and progression,
yet sex is rarely treated as a biological variable in hepatic immunology. This study
asks whether the liver immune composition itself differs by sex, and whether any
difference is driven by disease or is already present in disease-free liver. It is
built entirely on public data: public data in, published signatures out.

![Language](https://img.shields.io/badge/language-R-276DC3)
![License](https://img.shields.io/badge/license-MIT-green)

## Overview

The pipeline assigns each sample's sex directly from sex-chromosome gene expression
(so it works even where the metadata omits sex), scores immune cell-type and
functional-state signatures per sample, and estimates the sex effect for every
readout jointly across six public MASLD cohorts, adjusting for fibrosis stage. Every
headline result must survive a deliberately demanding gauntlet: two independent
fibrosis codings, false-discovery-rate correction, leave-one-cohort-out resampling,
a random-effects meta-analysis (with between-cohort heterogeneity I², and a
prediction interval), a lineage-specificity control (adjusting each subset for its
parent compartment), an independent deconvolution method, a metabolic (BMI /
type-2 diabetes / age) deconfounding step, and a disease-free (GTEx) control. The
resulting sex directions are then examined at single-cell resolution in an in-house
CITE-seq cohort and four independent public single-cell cohorts.

Two reproducible axes emerge. The first is a **baseline, male-biased MAIT-cell
signature** — already present in disease-free liver, specific to the receptor-defined
(*SLC4A10* / *TRAV1-2*) MAIT programme, homogeneous across cohorts (I² = 0%), and
*strengthening* rather than weakening once total T-cell content is accounted for. The
second is a **disease-associated, female-biased programme** — anchored on regulatory
T cells, with CD8 T cells and conventional dendritic cells — that is flat in
disease-free liver and appears only with MASLD.

## Repository structure

```
MASLD_sex_immune/
├── README.md
├── LICENSE
├── results_R/                 Regenerated outputs (figures and tables; git-ignored)
└── scripts/
    ├── run_all.R                  end-to-end runner (analysis, meta, tables, figures)
    ├── install_packages.R         install all dependencies
    ├── 00_config.R                shared paths, signatures and helpers (edit paths here)
    ├── geo_loaders.R              GEO download and parsing helpers
    ├── 01_sex_assignment.R        expression-based sex calls (XIST and Y-panel)
    ├── 02_build_matrix.R          singscore signature matrix and covariates
    ├── 03_maineffect.R            sex main-effect mixed models (FDR, LOCO)
    ├── 04_gtex_control.R          disease-free GTEx liver control
    ├── 05_export_expr.R           linear-scale matrices for deconvolution
    ├── 06_prep_deconv.R           symbol-level matrices for immunedeconv
    ├── 07_run_deconv.R            xCell and MCP-counter deconvolution
    ├── 08_concordance.R           deconvolution vs signature concordance
    ├── 09_deconfound.R            BMI / type-2 diabetes / age deconfounding
    ├── 13_meta_random.R           random-effects meta-analysis (metafor; Supp Fig S2/S8)
    ├── 15_lineage_specificity.R   subset vs parent-lineage adjustment
    ├── sc_build_andrews.R         build a Seurat object (Andrews 2024)
    ├── sc_build_guilliams.R       build a Seurat object (Guilliams / Liver Cell Atlas)
    ├── sc_build_ramachandran.R    build a Seurat object (Ramachandran 2019)
    ├── sc_validate_andrews.R      gate cell types and test sex directions (Andrews)
    ├── sc_validate_hlica.R        gate cell types and test sex directions (HLiCA)
    ├── sc_validate_ramachandran.R gate cell types and test sex directions (Ramachandran)
    ├── sc_guilliams_clean.R       patient-level re-analysis of the Guilliams cohort
    ├── sc_build_owncohort.R       in-house CITE-seq per-patient fractions
    ├── sc_mait_functional.R       male vs female MAIT functional-state comparison
    ├── sc_meta_forest.R           Cliff's delta single-cell effect sizes per cohort
    ├── 16_sc_pooled_mait.R        pooled per-donor single-cell test (mixed model + van Elteren)
    ├── make_tables.R              Tables 1-2 + Supplementary Tables S1-S5
    ├── fig1_sex_assignment.R      Figure 1 (sex assignment)
    ├── fig2_main_effects.R        Figure 2 (forest + distributions)
    ├── fig3_mait_identity.R       Figure 3 (MAIT receptor identity)
    ├── fig4_gtex_specificity.R    Figure 4 (MASLD vs GTEx + robustness)
    ├── fig5_deconvolution.R       Figure 5 (deconvolution concordance)
    ├── fig6_singlecell.R          Figure 6 (single-cell)
    ├── make_workflow_schematic.R  Supplementary Figure S1 (analysis workflow)
    ├── suppfig_deconfounding.R    Supplementary Figure S4 (metabolic deconfounding)
    ├── 11_supp_figures.R          Supplementary Figures S2, S3, S5-S8
    └── FIGURES_LIST.md            figure and table catalogue
```

## Requirements

R (>= 4.1) with the packages listed in `scripts/install_packages.R`, including
`data.table`, `lme4`, `lmerTest`, `metafor`, `singscore`, `immunedeconv` (xCell,
MCP-counter), `Seurat`, `ggplot2`, `fgsea`, `msigdbr` and `openxlsx`. Install them once with:

```r
source("scripts/install_packages.R")
```

For exact reproducibility, capture package versions with `sessionInfo()` after
installation, or manage them with [renv](https://rstudio.github.io/renv/)
(`renv::init()` then `renv::snapshot()`).

## Quick start

Run everything in dependency order (analysis, meta, tables, figures) from the
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
| `RUN_META`    | on      | random-effects meta-analysis (`13`) + lineage specificity (`15`) |
| `RUN_SC`      | off     | single-cell validations + pooled per-donor test (`16`) (require the large raw downloads) |
| `RUN_TABLES`  | on      | `make_tables.R` |
| `RUN_FIGURES` | on      | `fig1`–`fig6`, `make_workflow_schematic.R`, `suppfig_deconfounding.R`, `11_supp_figures.R` |

### Reproduce the bulk results without downloading raw data

The full pipeline regenerates `analysis_matrix.csv`, the deconvolution scores and the
per-donor single-cell fractions into the analysis output directory (`MASLD_REALDIR`).
Once those exist (either from a full run, or copied into `data/` for a release), the
downstream steps reproduce the core results directly:

```r
setwd("scripts")
Sys.setenv(MASLD_REALDIR = normalizePath("../results_R"))
source("03_maineffect.R"); source("13_meta_random.R"); source("15_lineage_specificity.R")
source("08_concordance.R"); source("09_deconfound.R"); source("make_tables.R")
```

### Full pipeline from raw data

Edit the paths in `00_config.R`, then run `run_all.R` (optionally with `RUN_SC=1` and
`RUN_INSTALL=1`), or run the numbered scripts `01`–`09`, `13`, `15`, the `sc_*`
scripts, `16`, `make_tables.R` and the figure scripts in that order.
`FIGURES_LIST.md` maps every figure and table to the script that produces it.

## Data availability

All primary data are public. Raw data are not redistributed here; the scripts
retrieve them from source, and `data/` holds only the small processed inputs needed
to reproduce the results. No controlled-access or patient-identifiable data are used;
the in-house CITE-seq file contains aggregate per-patient counts only. Before
archiving a versioned release, refresh `data/analysis_matrix.csv`, the deconvolution
score tables and the per-donor single-cell fraction files from the analysis output.

| Source | Accession |
|--------|-----------|
| MASLD bulk cohorts | GEO GSE130970, GSE162694, GSE135251, GSE89632, GSE167523, GSE48452 |
| Disease-free liver | GTEx v8 |
| Single-cell, Andrews 2024 | GEO GSE243981 |
| Single-cell, Ramachandran 2019 | GEO GSE136103 |
| Single-cell, Guilliams 2022 (Liver Cell Atlas) | GEO GSE192742 |
| Single-cell, HLiCA | CZ CELLxGENE (cellxgene.cziscience.com) |

## Citation

If you use this code, please cite the associated manuscript (in preparation). A
versioned release will be archived with a Zenodo DOI on acceptance. See
[CITATION.cff](CITATION.cff).

## License

Released under the MIT License. See [LICENSE](LICENSE).
