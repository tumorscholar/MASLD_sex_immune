## 16_sc_pooled_mait.R --------------------------------------------------------
## Better-powered single-cell test of the sex effect. sc_meta_forest.R meta-
## analyses four small, sex-imbalanced cohorts as separate Cliff's deltas, which
## for MAIT is underpowered and heterogeneous (pooled I2~85%, p~0.97). This script
## instead pools ALL per-donor fractions into a SINGLE analysis with cohort as a
## stratum / random intercept, which is far more powerful and is the standard way
## to combine per-subject data across small studies. Two complementary tests per
## readout (MAIT / Treg / CD8), so the conclusion does not hinge on distributional
## assumptions:
##   (1) linear mixed model:  asin(sqrt(fraction)) ~ sexM + (1 | cohort)
##   (2) van Elteren test:    cohort-stratified Wilcoxon (nonparametric backstop)
## Also reports per-cohort direction so consistency is explicit.
##
## Uses the SAME cohort registry / filters / columns as sc_meta_forest.R, so if
## that script runs, this one does too. Reads the per-donor fraction CSVs written
## by the single-cell validations (needs those outputs present).
## ===========================================================================
source("00_config.R")
suppressMessages({library(lme4); library(lmerTest)})
OUTBASE <- file.path(REALDIR, "single_cell")
MIN_T   <- 50

## cohort registry (identical to sc_meta_forest.R)
CO <- list(
 list(name="Andrews (healthy)",        csv="andrews2024/out/Andrews2024_per_donor_fractions.csv",      sub=c(group="healthy")),
 list(name="HLiCA (healthy)",          csv="hlica/out/HLiCA_per_patient.csv",                           sub=NULL),
 list(name="Guilliams (lean)",         csv="guilliams/out/Guilliams_clean_per_patient.csv",            sub=c(diet="Lean")),
 list(name="Ramachandran (cirrhotic)", csv="ramachandran/out/Ramachandran2019_per_donor_fractions.csv",sub=c(group="cirrhotic")),
 list(name="Ramachandran (healthy)",   csv="ramachandran/out/Ramachandran2019_per_donor_fractions.csv",sub=c(group="healthy"))
)
READOUTS <- c(MAIT_pctT="MAIT", Treg_pctT="Treg", CD8_pctT="CD8")

## ---- stack all per-donor rows into one long table ---------------------------
donors <- list()
for (co in CO) {
  fp <- file.path(OUTBASE, co$csv)
  if (!file.exists(fp)) { cat("missing:", co$csv, "\n"); next }
  x <- read.csv(fp, stringsAsFactors = FALSE)
  names(x)[names(x) == "n_Tcells"] <- "n_T"
  if (!is.null(co$sub)) x <- x[as.character(x[[names(co$sub)]]) == unname(co$sub), , drop = FALSE]
  if ("n_T" %in% names(x)) x <- x[!is.na(x$n_T) & x$n_T >= MIN_T, ]
  x <- x[x$sex %in% c("M","F"), ]
  if (!nrow(x)) next
  keep <- c("sex", intersect(names(READOUTS), names(x)))
  d <- x[, keep, drop = FALSE]
  d$n_T <- if ("n_T" %in% names(x)) x$n_T else NA_real_   # keep T-cell capture for the power check
  d$cohort <- co$name
  donors[[co$name]] <- d
}
D <- do.call(rbind, donors)
if (is.null(D) || !nrow(D)) stop("no single-cell per-donor CSVs found under ", OUTBASE)
cat(sprintf("pooled donors: %d across %d cohorts | M/F: %d/%d\n",
    nrow(D), length(unique(D$cohort)), sum(D$sex=="M"), sum(D$sex=="F")))

## ---- van Elteren cohort-stratified Wilcoxon (Z>0 = higher in males) ----------
van_elteren <- function(value, sex, stratum) {
  ok <- !is.na(value) & sex %in% c("M","F")
  value <- value[ok]; sex <- sex[ok]; stratum <- stratum[ok]
  num <- 0; den <- 0
  for (h in unique(stratum)) {
    i <- stratum == h; v <- value[i]; s <- sex[i]
    nh <- length(v); nM <- sum(s=="M"); nF <- sum(s=="F")
    if (nM < 1 || nF < 1) next
    w <- 1/(nh + 1); r <- rank(v)                       # mid-ranks within stratum
    num <- num + w * (sum(r[s=="M"]) - nM*(nh+1)/2)
    den <- den + w^2 * (nM*nF*(nh+1)/12)
  }
  if (den <= 0) return(c(Z=NA, p=NA))
  Z <- num/sqrt(den); c(Z=Z, p=2*pnorm(-abs(Z)))
}

## ---- per readout: pooled mixed model + van Elteren + per-cohort direction ----
out <- list()
for (rv in names(READOUTS)) {
  if (!rv %in% names(D)) next
  d <- D[!is.na(D[[rv]]), c("sex", rv, "cohort")]; names(d)[2] <- "y"
  d$sexM <- as.numeric(d$sex == "M")
  d$yt   <- asin(sqrt(pmin(pmax(d$y/100, 0), 1)))       # arcsine-sqrt of proportion
  ## (1) mixed model
  mm <- tryCatch(lmerTest::lmer(yt ~ sexM + (1|cohort), d,
        control = lmerControl(calc.derivs = FALSE)), error = function(e) NULL)
  if (is.null(mm)) {                                    # fall back to cohort fixed effect
    mm2 <- lm(yt ~ sexM + cohort, d); co <- summary(mm2)$coefficients
    b <- co["sexM","Estimate"]; p <- co["sexM","Pr(>|t|)"]; model <- "lm(+cohort FE)"
  } else { co <- summary(mm)$coefficients
    b <- co["sexM","Estimate"]; p <- co["sexM","Pr(>|t|)"]; model <- "lmer(1|cohort)" }
  ## (2) van Elteren
  ve <- van_elteren(d$y, d$sex, d$cohort)
  ## per-cohort direction (sign of median male - median female)
  dirs <- sapply(split(d, d$cohort), function(z) {
    if (sum(z$sexM==1) < 1 || sum(z$sexM==0) < 1) return(NA)
    sign(median(z$y[z$sexM==1]) - median(z$y[z$sexM==0])) })
  nmale_dir <- sum(dirs > 0, na.rm=TRUE); nfem_dir <- sum(dirs < 0, na.rm=TRUE)
  out[[rv]] <- data.frame(readout=READOUTS[[rv]], model=model,
    nM=sum(d$sexM==1), nF=sum(d$sexM==0),
    beta_arcsin=round(b,4), p_mixed=signif(p,3),
    vanElteren_Z=round(ve["Z"],3), p_vanElteren=signif(ve["p"],3),
    cohorts_male_dir=nmale_dir, cohorts_female_dir=nfem_dir, row.names=NULL)
}
res <- do.call(rbind, out)
write.csv(res, file.path(OUTBASE, "sc_pooled_perdonor.csv"), row.names=FALSE)

## ---- rare-cell power check: how many MAIT cells actually underlie each fraction? --
if ("n_T" %in% names(D) && "MAIT_pctT" %in% names(D)) {
  D$mait_cells <- D$n_T * D$MAIT_pctT / 100
  m <- D$mait_cells[!is.na(D$mait_cells)]
  cat("\n==== per-donor MAIT cell counts (rare-cell power check) ====\n")
  if (length(m)) cat(sprintf(
    "pooled (%d donors with n_T): median %.1f MAIT cells/donor (IQR %.1f-%.1f); %.0f%% have <5, %.0f%% have <10\n",
    length(m), median(m), quantile(m,.25), quantile(m,.75),
    100*mean(m < 5), 100*mean(m < 10)))
  bycoh <- do.call(rbind, lapply(split(D, D$cohort), function(z) {
    mm <- z$mait_cells[!is.na(z$mait_cells)]
    data.frame(cohort=z$cohort[1], donors=nrow(z),
      median_nT=if(all(is.na(z$n_T))) NA else round(median(z$n_T,na.rm=TRUE),0),
      median_MAIT_cells=if(!length(mm)) NA else round(median(mm),1),
      pct_under5=if(!length(mm)) NA else round(100*mean(mm<5),0), row.names=NULL) }))
  print(bycoh, row.names=FALSE)
  cat("If most donors carry only a handful of MAIT cells, per-donor MAIT fractions are\n",
      "noise-dominated and the pooled null is a POWER ceiling, not absence of effect -\n",
      "which is why the whole-tissue bulk signature is the more sensitive measure here.\n")
} else cat("\n(no n_T column available - cannot compute per-donor MAIT cell counts)\n")
cat("\n==== pooled per-donor single-cell sex effect (beta/Z > 0 = higher in MALES) ====\n")
print(res, row.names=FALSE)
cat("\nInterpretation: this pooled per-donor test uses every donor jointly with cohort\n",
    "as a stratum, rather than meta-analysing 4 tiny cohorts. If MAIT is male-biased\n",
    "here (p_mixed and/or p_vanElteren < 0.05) with most cohorts pointing the same way,\n",
    "it is a stronger single-cell confirmation than the heterogeneous cohort-level pool.\n")
cat("Wrote sc_pooled_perdonor.csv to", OUTBASE, "\n")
