## 15_lineage_specificity.R ---------------------------------------------------
## Specificity control for the sex-divergent immune readouts. A reviewer can
## argue the female-biased programme is just "females have more immune infiltrate"
## and the male-MAIT bias just "more T cells". This script answers that directly:
## for each headline subset it re-fits the sex effect ADJUSTED for its lineage
## denominator (T subsets vs total T cells; myeloid subsets vs monocyte/macrophage;
## B subsets vs total B), and reports how much the sex coefficient attenuates.
##
##   fit0:  readout ~ sexM + fibrosis + (1|cohort)              (unadjusted)
##   fit1:  readout ~ sexM + fibrosis + lineage + (1|cohort)    (lineage-adjusted)
##
## A subset that stays significant and same-sign after adjustment is a COMPOSITIONAL
## shift specific to that subset, not explained by total lineage abundance. Note the
## key asymmetry: total T cells (ct_Tcell) are themselves female-biased, so adjusting
## the MALE-biased MAIT readouts for total T is a hard test (it removes a female-
## leaning denominator) - surviving it is strong evidence the male-MAIT bias is real.
## Self-contained: reads analysis_matrix.csv only.
## ===========================================================================
source("00_config.R")
suppressMessages({library(lme4); library(lmerTest)})

mtx <- read.csv(file.path(REALDIR, "analysis_matrix.csv"), stringsAsFactors = FALSE)
mtx$fib_ord <- code_fibrosis(mtx$fibrosis_stage)
mtx <- mtx[mtx$sex_assigned %in% c("F","M") & !is.na(mtx$fib_ord), , drop = FALSE]
mtx$sexM <- as.numeric(mtx$sex_assigned == "M")
cat(sprintf("samples: %d | cohorts: %d | F/M: %d/%d\n",
    nrow(mtx), length(unique(mtx$cohort)), sum(mtx$sexM==0), sum(mtx$sexM==1)))

## within-cohort z-score every readout + the lineage denominators (comparable SD units)
READOUTS <- grep("^(ct_|st_)", names(mtx), value = TRUE)
for (r in READOUTS) mtx[[r]] <- ave(mtx[[r]], mtx$cohort, FUN = zscore)

## lineage denominator for each subset (its parent compartment)
Tsub  <- c("ct_MAIT","ct_MAITspec","ct_MAITpromisc","ct_Treg","ct_CD8T",
           "st_Th1","st_Th17","st_exhaustion","st_Tpex","st_Trm","st_cytotoxCD4","st_CD8_cytotox")
Mysub <- c("ct_DC","ct_Neutrophil")
Bsub  <- c("ct_Plasma")
denom_of <- function(r) if (r %in% Tsub) "ct_Tcell" else
                        if (r %in% Mysub) "ct_MonoMac" else
                        if (r %in% Bsub) "ct_Bcell" else NA_character_

beta_sex <- function(form) {
  m <- tryCatch(lmerTest::lmer(form, mtx, REML = TRUE,
        control = lmerControl(calc.derivs = FALSE)), error = function(e) NULL)
  if (is.null(m)) return(c(beta = NA, p = NA))
  co <- summary(m)$coefficients
  if (!"sexM" %in% rownames(co)) return(c(beta = NA, p = NA))
  c(beta = co["sexM","Estimate"], p = co["sexM","Pr(>|t|)"])
}

## context: are the lineage denominators themselves sex-biased?
cat("\n== lineage denominators: own sex effect (beta<0 = higher in females) ==\n")
for (d in c("ct_Tcell","ct_MonoMac","ct_Bcell")) {
  b <- beta_sex(as.formula(sprintf("`%s` ~ sexM + fib_ord + (1|cohort)", d)))
  cat(sprintf("  %-12s beta=%+.3f  p=%.2g  (%s)\n", d, b["beta"], b["p"],
              ifelse(b["beta"] > 0, "male", "female")))
}

HEAD <- c("ct_MAIT","ct_MAITspec","ct_MAITpromisc","ct_Treg","ct_CD8T","st_Th1",
          "st_exhaustion","st_Tpex","ct_DC","ct_Neutrophil","ct_Plasma")
rows <- list()
for (r in HEAD) {
  if (!r %in% names(mtx)) next
  den <- denom_of(r)
  b0 <- beta_sex(as.formula(sprintf("`%s` ~ sexM + fib_ord + (1|cohort)", r)))
  if (is.na(den)) { b1 <- c(beta = NA, p = NA) } else
    b1 <- beta_sex(as.formula(sprintf("`%s` ~ sexM + fib_ord + `%s` + (1|cohort)", r, den)))
  atten <- if (is.na(b1["beta"]) || b0["beta"] == 0) NA else round(100 * (1 - b1["beta"]/b0["beta"]), 1)
  rows[[r]] <- data.frame(readout = r, denom = den,
    beta_unadj = round(b0["beta"],3), p_unadj = signif(b0["p"],3),
    beta_adj = round(b1["beta"],3),  p_adj = signif(b1["p"],3),
    pct_atten = atten, survives = !is.na(b1["p"]) & b1["p"] < 0.05 &
      sign(b1["beta"]) == sign(b0["beta"]), row.names = NULL)
}
res <- do.call(rbind, rows)
res$fdr_adj <- signif(p.adjust(res$p_adj, method = "BH"), 3)
write.csv(res, file.path(REALDIR, "lineage_specificity.csv"), row.names = FALSE)

cat("\n==== sex effect BEFORE vs AFTER adjusting for lineage denominator ====\n")
cat("(survives = adjusted p<0.05 AND same direction as unadjusted)\n")
print(res, row.names = FALSE)

cat("\nInterpretation:\n")
cat(" - MAIT / MAITspec are MALE-biased; total T cells are FEMALE-biased, so",
    "\n   adjusting MAIT for total T is a stringent test. If the male effect survives,",
    "\n   the male-MAIT bias is a within-T compositional shift, not more T cells.\n")
cat(" - Treg / CD8 / Th1 staying female-biased after adjusting for total T means",
    "\n   the female programme is subset-specific, not generic immune-richness.\n")
cat("Wrote lineage_specificity.csv to", REALDIR, "\n")
