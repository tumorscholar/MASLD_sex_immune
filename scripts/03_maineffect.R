## 03_maineffect.R -----------------------------------------------------------
## Sex MAIN-EFFECT scan (stage-adjusted), the core test of the paper.
##   readout ~ sexM + fibrosis + (1 | cohort)
## Fits the sex effect for every immune readout, adjusting for fibrosis stage
## and treating cohort as a random intercept (internal meta-analysis).
## BH-FDR across readouts; leave-one-cohort-out; both fibrosis codings;
## metabolic-adjusted sub-analysis.
##
## Fixed-effect p-values use lmerTest (Satterthwaite degrees of freedom).
## ---------------------------------------------------------------------------
source("00_config.R")
suppressMessages({library(lme4); library(lmerTest)})

df <- read.csv(file.path(REALDIR, "analysis_matrix.csv"), stringsAsFactors = FALSE)
READOUTS <- grep("^(ct_|st_)", names(df), value = TRUE)

df$fib_ord <- code_fibrosis(df$fibrosis_stage)
df$fib_bin <- as.numeric(df$fib_ord > 0)
df <- df[df$sex_assigned %in% c("F","M"), , drop = FALSE]
df$sexM <- as.numeric(df$sex_assigned == "M")
df <- df[!is.na(df$fib_ord), , drop = FALSE]

## z-score each readout WITHIN cohort (comparable betas across studies)
for (r in READOUTS)
  df[[r]] <- ave(df[[r]], df$cohort, FUN = function(x) zscore(x))

cat("Samples after QC:", nrow(df),
    "| cohorts:", paste(names(table(df$cohort)), table(df$cohort), sep="=", collapse=" "), "\n")
cat("Sex counts (F/M):", sum(df$sexM == 0), "/", sum(df$sexM == 1), "\n")

## fit sexM effect from a mixed model; returns c(beta, se, p) or NAs
fit_main <- function(d, r, coding) {
  d <- d[stats::complete.cases(d[, c(r, coding, "sexM")]), , drop = FALSE]
  f <- as.formula(sprintf("%s ~ sexM + %s + (1|cohort)", r, coding))
  out <- tryCatch({
    m  <- lmerTest::lmer(f, data = d, REML = TRUE,
                         control = lmerControl(calc.derivs = FALSE))
    co <- summary(m)$coefficients
    c(co["sexM","Estimate"], co["sexM","Std. Error"], co["sexM","Pr(>|t|)"])
  }, error = function(e) c(NA, NA, NA))
  out
}

rows <- list()
for (coding in c("fib_ord","fib_bin")) {
  for (r in READOUTS) {
    d <- df[stats::complete.cases(df[, c(r, coding, "sexM")]), , drop = FALSE]
    est <- fit_main(d, r, coding); b <- est[1]; se <- est[2]; p <- est[3]
    ## leave-one-cohort-out
    signs <- c(); nsig <- 0
    for (co in unique(d$cohort)) {
      e2 <- fit_main(d[d$cohort != co, , drop = FALSE], r, coding)
      if (!is.na(e2[1])) { signs <- c(signs, sign(e2[1])); nsig <- nsig + (e2[3] < 0.05) }
    }
    loco_consistent <- length(signs) > 0 && !is.na(b) && all(signs == sign(b))
    rows[[length(rows)+1]] <- data.frame(
      readout = r, coding = coding, n = nrow(d),
      beta_sexM = b, se = se, p = p,
      loco_sign_consistent = loco_consistent, loco_n_sig = nsig,
      stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)
res$fdr <- NA_real_
ok <- !is.na(res$p)
res$fdr[ok] <- p.adjust(res$p[ok], method = "BH")

## concordance across the two fibrosis codings
conc <- list()
for (r in READOUTS) {
  bo <- res$beta_sexM[res$readout==r & res$coding=="fib_ord"]
  bb <- res$beta_sexM[res$readout==r & res$coding=="fib_bin"]
  fo <- res$fdr[res$readout==r & res$coding=="fib_ord"]
  fb <- res$fdr[res$readout==r & res$coding=="fib_bin"]
  if (length(bo) && length(bb) && !is.na(bo) && !is.na(bb)) {
    conc[[length(conc)+1]] <- data.frame(
      readout = r, beta_ord = bo, beta_bin = bb, fdr_ord = fo, fdr_bin = fb,
      same_sign = sign(bo) == sign(bb),
      headline = (!is.na(fo) && !is.na(fb) && fo < 0.10 && fb < 0.10 && sign(bo) == sign(bb)),
      stringsAsFactors = FALSE)
  }
}
conc <- do.call(rbind, conc)
conc <- conc[order(conc$fdr_ord), ]

write.csv(res[order(res$fdr), ], file.path(REALDIR, "maineffect_results.csv"), row.names = FALSE)
write.csv(conc, file.path(REALDIR, "maineffect_concordance.csv"), row.names = FALSE)

## metabolic-adjusted sex main effect (samples with BMI + T2D) – plain OLS
g <- df[stats::complete.cases(df[, c("bmi","t2d")]), , drop = FALSE]
if (nrow(g) > 20) {
  cat("\nMetabolic-adjusted sex main effect on n=", nrow(g),
      " (", paste(unique(g$cohort), collapse=","), ")\n", sep="")
  mrows <- list()
  for (r in READOUTS) {
    d <- g[!is.na(g[[r]]), , drop = FALSE]
    fit <- tryCatch(lm(as.formula(sprintf("%s ~ sexM + fib_ord + bmi + t2d + age", r)), d),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      s <- summary(fit)$coefficients
      mrows[[length(mrows)+1]] <- data.frame(readout = r,
        beta_sexM = s["sexM","Estimate"], p_sexM = s["sexM","Pr(>|t|)"], stringsAsFactors = FALSE)
    }
  }
  write.csv(do.call(rbind, mrows), file.path(REALDIR, "maineffect_metabolic.csv"), row.names = FALSE)
}

cat("\n===== TOP sex MAIN effects (stage-adjusted, by FDR) =====\n")
print(head(res[order(res$fdr),
    c("readout","coding","n","beta_sexM","p","fdr","loco_sign_consistent","loco_n_sig")], 12),
    row.names = FALSE)
cat("\n===== CONCORDANCE across codings (headline = FDR<0.10 both + same sign) =====\n")
print(conc, row.names = FALSE)
head_hits <- conc$readout[conc$headline]
cat("\nHEADLINE sex-different readouts:",
    if (length(head_hits)) paste(head_hits, collapse=", ") else "none", "\n")
cat("(beta_sexM > 0 = higher in MALES; < 0 = higher in FEMALES)\n")
cat("Wrote: maineffect_results.csv, maineffect_concordance.csv, maineffect_metabolic.csv\n")
