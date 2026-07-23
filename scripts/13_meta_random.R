## 13_meta_random.R -----------------------------------------------------------
## Formal random-effects meta-analysis of the sex effect, for both the bulk
## cohorts and the single-cell cohorts. The core analysis (03_maineffect.R) pools
## cohorts with a random intercept; this script complements it with an explicit
## per-cohort effect-size meta-analysis (metafor, REML), reporting between-cohort
## heterogeneity (I-squared, tau-squared, Q) and 95% prediction intervals, and
## draws standard forest plots. This is the presentation a meta-analysis reviewer
## expects and makes the heterogeneity of each effect explicit.
##
## Bulk: per-cohort sex coefficient (readout ~ sexM + fibrosis) on within-cohort
##       z-scored readouts, pooled by metafor::rma.
## Single-cell: per-cohort Cliff's delta + SE from sc_meta_forest.R, pooled the
##       same way. Note: Ramachandran contributes disjoint healthy and cirrhotic
##       subsets as separate rows; they are flagged in the output.
## ===========================================================================
source("00_config.R")
suppressMessages(library(metafor))
FIG <- FIGDIR; dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
MALE <- "#2c7fb8"; FEM <- "#d95f0e"
save_fig <- function(draw, name, w, h) {
  pdf (file.path(FIG, paste0(name,".pdf")),  width=w, height=h);                                        draw(); dev.off()
  tiff(file.path(FIG, paste0(name,".tiff")), width=w, height=h, units="in", res=300, compression="lzw"); draw(); dev.off()
  png (file.path(FIG, paste0(name,".png")),  width=w, height=h, units="in", res=300);                   draw(); dev.off()
  cat("wrote", name, "(pdf, tiff, png)\n")
}
LAB <- c(ct_MAIT="MAIT", ct_MAITspec="MAIT (SLC4A10/TRAV1-2)", ct_Treg="Treg",
         ct_CD8T="CD8 T", ct_DC="Dendritic (cDC)", st_Th1="Th1")
HEAD <- names(LAB)

## ===================== BULK random-effects meta =====================
mtx <- tryCatch(read.csv(file.path(REALDIR, "analysis_matrix.csv"), stringsAsFactors=FALSE),
                error=function(e) NULL)
if (!is.null(mtx)) {
  READOUTS <- grep("^(ct_|st_)", names(mtx), value=TRUE)
  mtx$fib_ord <- code_fibrosis(mtx$fibrosis_stage)
  mtx$sexM    <- as.numeric(mtx$sex_assigned == "M")
  mtx <- mtx[mtx$sex_assigned %in% c("F","M") & !is.na(mtx$fib_ord), , drop=FALSE]
  ## within-cohort z-score so per-cohort betas are in comparable SD units
  for (r in READOUTS) mtx[[r]] <- ave(mtx[[r]], mtx$cohort, FUN=function(x) zscore(x))

  per_cohort <- function(r) {
    do.call(rbind, lapply(unique(mtx$cohort), function(g) {
      d <- mtx[mtx$cohort==g, ]; d <- d[stats::complete.cases(d[, c(r,"sexM","fib_ord")]), ]
      if (nrow(d) < 8 || length(unique(d$sexM)) < 2) return(NULL)
      fit <- tryCatch(lm(as.formula(sprintf("`%s` ~ sexM + fib_ord", r)), d), error=function(e) NULL)
      if (is.null(fit)) return(NULL)
      co <- summary(fit)$coefficients
      if (!"sexM" %in% rownames(co)) return(NULL)
      data.frame(cohort=g, yi=co["sexM","Estimate"], sei=co["sexM","Std. Error"],
                 n=nrow(d), row.names=NULL)
    }))
  }
  brows <- list()
  for (r in READOUTS) {
    pc <- per_cohort(r); if (is.null(pc) || nrow(pc) < 2) next
    m <- tryCatch(rma(yi=pc$yi, sei=pc$sei, method="REML"), error=function(e) NULL)
    if (is.null(m)) next
    pr <- predict(m)
    brows[[r]] <- data.frame(readout=r, k=m$k, beta=as.numeric(m$b), se=m$se,
      ci_lo=m$ci.lb, ci_hi=m$ci.ub, pval=m$pval,
      I2=m$I2, tau2=m$tau2, Q=m$QE, Q_p=m$QEp,
      PI_lo=pr$pi.lb, PI_hi=pr$pi.ub, row.names=NULL)
  }
  bmeta <- do.call(rbind, brows)
  if (!is.null(bmeta)) { bmeta <- bmeta[order(bmeta$pval), ]
    write.csv(bmeta, file.path(REALDIR, "bulk_meta_random.csv"), row.names=FALSE)
    cat("\n==== BULK random-effects meta (beta>0 = higher in males) ====\n")
    print(bmeta[, c("readout","k","beta","ci_lo","ci_hi","pval","I2","Q_p","PI_lo","PI_hi")],
          row.names=FALSE, digits=3)
  }
  ## forest plots for the headline readouts
  hd <- HEAD[HEAD %in% READOUTS]
  if (length(hd)) save_fig(function() {
    par(mfrow=c(length(hd),1), mar=c(4,4,2,2))
    for (r in hd) { pc <- per_cohort(r); if (is.null(pc) || nrow(pc) < 2) { plot.new(); next }
      m <- rma(yi=pc$yi, sei=pc$sei, method="REML")
      forest(m, slab=pc$cohort, xlab="sex effect (SD; >0 male)", header=LAB[[r]],
             addpred=TRUE, cex=0.8, col=MALE, mlab=sprintf("RE model (I2=%.0f%%)", m$I2))
    }
  }, "SuppFig_S5_bulk_meta_forest", 7, 2.2*length(hd))
}

## ===================== SINGLE-CELL random-effects meta =====================
SC <- file.path(REALDIR, "single_cell")
eff <- tryCatch(read.csv(file.path(SC, "meta/sc_meta_effects.csv"), stringsAsFactors=FALSE),
                error=function(e) NULL)
if (!is.null(eff)) {
  eff <- eff[eff$cohort != "POOLED" & is.finite(eff$delta) & is.finite(eff$se) & eff$se > 0, ]
  srows <- list()
  for (rd in unique(eff$readout)) {
    s <- eff[eff$readout==rd, ]; if (nrow(s) < 2) next
    m <- tryCatch(rma(yi=s$delta, sei=s$se, method="REML"), error=function(e) NULL); if (is.null(m)) next
    pr <- predict(m)
    srows[[rd]] <- data.frame(readout=rd, k=m$k, delta=as.numeric(m$b), se=m$se,
      ci_lo=m$ci.lb, ci_hi=m$ci.ub, pval=m$pval, I2=m$I2, tau2=m$tau2, Q_p=m$QEp,
      PI_lo=pr$pi.lb, PI_hi=pr$pi.ub, row.names=NULL)
  }
  smeta <- do.call(rbind, srows)
  if (!is.null(smeta)) {
    write.csv(smeta, file.path(REALDIR, "sc_meta_random.csv"), row.names=FALSE)
    cat("\n==== SINGLE-CELL random-effects meta (Cliff's delta, >0 = male) ====\n")
    print(smeta[, c("readout","k","delta","ci_lo","ci_hi","pval","I2","PI_lo","PI_hi")],
          row.names=FALSE, digits=3)
    save_fig(function() {
      rds <- unique(eff$readout); par(mfrow=c(length(rds),1), mar=c(4,4,2,2))
      for (rd in rds) { s <- eff[eff$readout==rd, ]; if (nrow(s) < 2) { plot.new(); next }
        m <- rma(yi=s$delta, sei=s$se, method="REML")
        forest(m, slab=s$cohort, xlab="Cliff's delta (>0 male)", header=rd, addpred=TRUE,
               cex=0.8, col=MALE, mlab=sprintf("RE model (I2=%.0f%%)", m$I2)) }
    }, "SuppFig_S6_sc_meta_forest_RE", 7, 2.2*length(unique(eff$readout)))
  }
}
cat("\nDone. Wrote bulk_meta_random.csv, sc_meta_random.csv and forest figures to",
    REALDIR, "\n")
