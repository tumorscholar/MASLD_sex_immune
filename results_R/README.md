# results_R

Regenerated analysis outputs. Most of this folder is git-ignored because it is
rebuilt by the scripts; only two things are tracked:

- `gtex_healthy_sex.csv` – the disease-free GTEx liver sex-effect table
  (deposited so Table 3 / Fig 5 reproduce without re-downloading GTEx v8).
- this README.

Running `scripts/run_all.R` (or the individual scripts) repopulates:

- `figures/`  – Fig 1-8 and Supp Fig S1-S3 as PDF / TIFF / PNG
- `tables/`   – Table 1-4 as CSV plus `MASLD_sex_tables.xlsx`
- the pipeline CSVs (`analysis_matrix.csv`, `maineffect_results.csv`, etc.)
