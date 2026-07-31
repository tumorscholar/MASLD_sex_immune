## 09_deconfound.R -----------------------------------------------------------
## Within-GSE89632 deconfounding test: does the sex effect survive BMI/T2D/age
## adjustment on the SAME samples? Fitting on one cohort isolates the effect of
## *adjustment* from the effect of *losing samples*, so any change in the sex
## beta can only mean metabolic confounding, not lost power.
## ---------------------------------------------------------------------------
source("00_config.R")

df <- read.csv(file.path(REALDIR, "analysis_matrix.csv"), stringsAsFactors = FALSE)
READOUTS <- grep("^(ct_|st_)", names(df), value = TRUE)
HEAD <- c("ct_MAIT","ct_Treg","ct_CD8T","ct_MonoMac","ct_DC","st_Th1","st_exhaustion")

df$fib_ord <- code_fibrosis(df$fibrosis_stage)
df <- df[df$cohort == "GSE89632" & df$sex_assigned %in% c("F","M"), , drop = FALSE]
df$sexM <- as.numeric(df$sex_assigned == "M")

rows <- list()
for (r in READOUTS) {
  d <- df[stats::complete.cases(df[, c(r,"fib_ord","bmi","t2d","age","sexM")]), , drop = FALSE]
  if (nrow(d) < 20) next
  d[[r]] <- zscore(d[[r]])   # z within these samples => comparable betas
  u <- tryCatch(lm(as.formula(sprintf("%s ~ sexM + fib_ord", r)), d), error = function(e) NULL)
  a <- tryCatch(lm(as.formula(sprintf("%s ~ sexM + fib_ord + bmi + t2d + age", r)), d), error = function(e) NULL)
  if (is.null(u) || is.null(a)) next
  su <- summary(u)$coefficients; sa <- summary(a)$coefficients
  bu <- su["sexM","Estimate"]; ba <- sa["sexM","Estimate"]
  rows[[length(rows)+1]] <- data.frame(
    readout = r, n = nrow(d), headline = r %in% HEAD,
    beta_unadj = round(bu,3), p_unadj = round(su["sexM","Pr(>|t|)"],4),
    beta_adj   = round(ba,3), p_adj   = round(sa["sexM","Pr(>|t|)"],4),
    pct_attenuation = if (bu != 0) round(100*(bu-ba)/bu,1) else NA_real_,
    stringsAsFactors = FALSE)
}
if (!length(rows)) {
  av <- colSums(!is.na(df[, c("bmi","t2d","age")]))
  stop(sprintf(paste0("09_deconfound: no GSE89632 samples have complete BMI/T2D/age, so the ",
    "deconfounding test cannot run. Non-NA counts among %d GSE89632 rows in analysis_matrix.csv - ",
    "bmi:%d  t2d:%d  age:%d. These covariates are not always exposed in GEO metadata; supply them ",
    "(e.g. merge a per-patient metabolic table into analysis_matrix.csv) to reproduce the ",
    "metabolic-deconfounding check (Results 3.4 / Supp Fig S4)."),
    nrow(df), av[["bmi"]], av[["t2d"]], av[["age"]]))
}
res <- do.call(rbind, rows)
res <- res[order(-res$headline, res$p_unadj), ]
write.csv(res, file.path(REALDIR, "deconf_within_cohort.csv"), row.names = FALSE)

cat("Within-GSE89632 (n<=58): sex effect UNADJUSTED vs BMI/T2D/age-ADJUSTED, same samples.\n")
cat("Small pct_attenuation + stable beta => metabolic state is NOT confounding the sex effect\n")
cat("(so any loss of significance vs the pooled n=616 result is power, not confounding).\n")
cat("Large attenuation (beta shrinks toward 0) => metabolic state explains it.\n\n")
print(res, row.names = FALSE)
cat("\n--- headline readouts only ---\n")
print(res[res$headline, ], row.names = FALSE)
cat("\nWrote deconf_within_cohort.csv\n")
