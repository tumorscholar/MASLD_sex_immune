## 08_concordance.R -----------------------------------------------------------
## Do the benchmarked deconvolution tools (xCell, MCP-counter) agree with the
## singscore sex DIRECTIONS? Re-fit the same mixed model on the deconvolution
## scores and compare the sign of the sex effect to the singscore reference.
## We test agreement in DIRECTION, not magnitude (methods never match on scale).
## ---------------------------------------------------------------------------
source("00_config.R")
suppressMessages({library(lme4); library(ggplot2)})
MALE_COL <- "#2c7fb8"; FEM_COL <- "#d95f0e"

meta <- read.csv(file.path(REALDIR, "analysis_matrix.csv"), stringsAsFactors = FALSE)
meta$fib_ord <- code_fibrosis(meta$fibrosis_stage)
meta$sexM <- ifelse(meta$sex_assigned == "M", 1,
              ifelse(meta$sex_assigned == "F", 0, NA))
meta <- meta[!is.na(meta$sexM), c("gsm","cohort","sexM","fib_ord")]

## singscore reference betas (ordinal coding)
me <- read.csv(file.path(REALDIR, "maineffect_results.csv"), stringsAsFactors = FALSE)
me <- me[me$coding == "fib_ord", ]
sing <- setNames(me$beta_sexM, me$readout)

## singscore readout -> deconvolution cell-type name maps
XMAP <- list(ct_Treg=c("Tregs"), ct_CD8T=c("CD8+ T-cells"), ct_DC=c("DC","cDC","pDC"),
             st_Th1=c("Th1 cells"), ct_NK=c("NK cells"),
             ct_MonoMac=c("Monocytes","Macrophages"), ct_Bcell=c("B-cells"),
             ct_Tcell=c("CD4+ T-cells"))
MMAP <- list(ct_CD8T=c("CD8 T cells"), ct_DC=c("Myeloid dendritic cells"),
             ct_NK=c("NK cells"), ct_MonoMac=c("Monocytic lineage"),
             ct_Bcell=c("B lineage"), ct_Tcell=c("T cells"),
             st_CD8_cytotox=c("Cytotoxic lymphocytes"))

## fit sexM sign on one deconvolution cell type
fit_ct <- function(df, col) {
  d <- df[stats::complete.cases(df[, c(col,"fib_ord","sexM")]), , drop = FALSE]
  if (length(unique(d[[col]])) < 5 || nrow(d) < 25) return(NULL)
  d$y <- zscore(d[[col]])
  m <- tryCatch(lme4::lmer(y ~ sexM + fib_ord + (1|cohort), data = d, REML = FALSE,
                           control = lmerControl(calc.derivs = FALSE)),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- lme4::fixef(m)["sexM"]
  list(beta = unname(b), n = nrow(d))
}

run <- function(fname, mapping, method) {
  fp <- file.path(REALDIR, fname)
  if (!file.exists(fp)) { message("missing ", fp); return(NULL) }
  dv <- read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE)
  if ("cohort" %in% names(dv)) dv$cohort <- NULL     # avoid cohort_x/cohort_y collision
  m <- merge(meta, dv, by.x = "gsm", by.y = "sample")
  cat(sprintf("%s: merged %d/%d samples; cell types=%d\n", method, nrow(m), nrow(meta), ncol(dv)-1))
  rows <- list()
  for (sig in names(mapping)) {
    sbeta <- if (sig %in% names(sing)) sing[[sig]] else NA
    sdir <- if (!is.na(sbeta)) (if (sbeta > 0) "MALE" else "female") else "-"
    for (ct in mapping[[sig]]) {
      if (!ct %in% names(m)) next
      r <- fit_ct(m, ct); if (is.null(r)) next
      ddir <- if (r$beta > 0) "MALE" else "female"
      agree <- if (is.na(sbeta)) "" else (if (ddir == sdir) "AGREE" else "DISAGREE")
      rows[[length(rows)+1]] <- data.frame(method = method, singscore_readout = sig,
        singscore_beta = if (!is.na(sbeta)) round(sbeta,3) else NA,
        singscore_dir = sdir, deconv_celltype = ct, deconv_beta = round(r$beta,3),
        n = r$n, deconv_dir = ddir, concordance = agree, stringsAsFactors = FALSE)
    }
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}

xc <- run("deconv_xcell.csv", XMAP, "xCell")
mc <- run("deconv_mcp.csv",   MMAP, "MCP-counter")
res <- rbind(xc, mc)
write.csv(res, file.path(REALDIR, "deconv_concordance.csv"), row.names = FALSE)
cat("\n==== DECONVOLUTION vs SINGSCORE sex directions ====\n")
print(res, row.names = FALSE)
prim <- res[res$concordance != "", ]
if (nrow(prim)) {
  ag <- sum(prim$concordance == "AGREE")
  cat(sprintf("\nDirection agreement: %d/%d matched cell types (%.0f%%)\n",
              ag, nrow(prim), 100*ag/nrow(prim)))
}

## ---- figure: singscore beta vs xCell beta ----
plotd <- xc[!is.na(xc$singscore_beta) & !is.na(xc$deconv_beta), ]
if (nrow(plotd)) {
  lim <- max(0.6, max(abs(c(plotd$singscore_beta, plotd$deconv_beta))) * 1.2)
  plotd$col <- ifelse(plotd$singscore_beta > 0, MALE_COL, FEM_COL)
  p <- ggplot(plotd, aes(singscore_beta, deconv_beta)) +
    annotate("rect", xmin=0, xmax=lim, ymin=0, ymax=lim, fill="#eef5ee") +
    annotate("rect", xmin=-lim, xmax=0, ymin=-lim, ymax=0, fill="#eef5ee") +
    geom_hline(yintercept=0, linewidth=.3) + geom_vline(xintercept=0, linewidth=.3) +
    geom_point(colour=plotd$col, size=3, shape=21, fill=plotd$col, stroke=.4) +
    geom_text(aes(label=deconv_celltype), size=2.5, hjust=-0.1, vjust=-0.4) +
    coord_cartesian(xlim=c(-lim,lim), ylim=c(-lim,lim)) +
    labs(x="singscore sex beta (>0 male)", y="xCell sex beta (>0 male)",
         title="Benchmarked deconvolution reproduces singscore sex directions") +
    theme_classic(base_size = 9)
  ggsave(file.path(FIGDIR, "Fig6_deconv_concordance.png"), p, width=5.6, height=5.2, dpi=300)
  cat("\nwrote figures/Fig6_deconv_concordance.png\n")
}
cat("\nNOTE: MAIT is not represented in xCell/MCP-counter references (TCR-defined, not in\n")
cat("bulk deconvolution signatures); its confirmation rests on the receptor-identity test,\n")
cat("the GTEx constitutional result, and the single-cell atlas. Wrote deconv_concordance.csv\n")
