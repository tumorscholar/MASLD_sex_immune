## 11_supp_figures.R ----------------------------------------------------------
## Candidate SUPPLEMENTARY figures, regenerated from the pipeline / single-cell
## output CSVs and saved as PDF (vector) + TIFF (300 dpi, LZW) + PNG.
## You decide later which of these become main vs supplementary.
##   S1  per-cohort single-cell validation dot plots (one per cohort)
##   S2  single-cell meta-analysis forest (MAIT / Treg / CD8)
##   S3  Guilliams QC – immune capture by sex + sex x diet
## ===========================================================================
RD  <- Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta")
SC  <- file.path(RD, "single_cell")
FIG <- file.path(RD, "figures"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
MALE <- "#2c7fb8"; FEM <- "#d95f0e"; MIN_T <- 50
rd <- function(p) if (file.exists(p)) read.csv(p, stringsAsFactors=FALSE) else { cat("missing:",p,"\n"); NULL }
save_fig <- function(draw, name, w, h) {
  render <- function() if (inherits(draw,"ggplot")) print(draw) else draw()
  pdf (file.path(FIG, paste0(name,".pdf")),  width=w, height=h);                                        render(); dev.off()
  tiff(file.path(FIG, paste0(name,".tiff")), width=w, height=h, units="in", res=300, compression="lzw"); render(); dev.off()
  png (file.path(FIG, paste0(name,".png")),  width=w, height=h, units="in", res=300);                   render(); dev.off()
  cat("wrote", name, "(pdf, tiff, png)\n")
}

## ---- S1: per-cohort validation dot plots -----------------------------------
COH <- list(
  list(id="Andrews",      csv="andrews2024/out/Andrews2024_per_donor_fractions.csv",       sub=c(group="healthy")),
  list(id="HLiCA",        csv="hlica/out/HLiCA_per_patient.csv",                            sub=NULL),
  list(id="Guilliams",    csv="guilliams/out/Guilliams_clean_per_patient.csv",             sub=c(diet="Lean")),
  list(id="Ramachandran", csv="ramachandran/out/Ramachandran2019_per_donor_fractions.csv", sub=NULL)
)
RMAP <- c(MAIT_pctT="MAIT (% T)", Treg_pctT="Treg (% T)", CD8_pctT="CD8 (% T)", cDC_pctImm="cDC (% imm)")
dotpanel <- function(x, cols, title) function() {
  par(mfrow=c(1,length(cols)), mar=c(3,3.4,3.5,0.5), oma=c(0,0,2,0), cex=0.8)
  isF <- x$sex=="F"; isM <- x$sex=="M"
  for (rv in cols) {
    f <- x[[rv]][isF]; m <- x[[rv]][isM]; p <- tryCatch(wilcox.test(m,f)$p.value, error=function(e) NA)
    plot(0,0,type="n", xlim=c(-0.6,1.6), ylim=range(c(f,m),na.rm=TRUE), xaxt="n", xlab="",
         ylab=if(rv==cols[1])"per-donor fraction (%)" else "", main=sprintf("%s\np=%.2g", RMAP[rv], p))
    points(jitter(rep(0,length(f)),amount=.09), f, pch=19, cex=.8, col=adjustcolor(FEM,.7))
    points(jitter(rep(1,length(m)),amount=.09), m, pch=19, cex=.8, col=adjustcolor(MALE,.7))
    segments(-0.25, median(f,na.rm=TRUE), 0.25, median(f,na.rm=TRUE), lwd=1.6)
    segments(0.75, median(m,na.rm=TRUE), 1.25, median(m,na.rm=TRUE), lwd=1.6)
    axis(1, at=c(0,1), labels=c(sprintf("F\n(%d)",length(f)), sprintf("M\n(%d)",length(m))), padj=0.5)
  }
  mtext(paste0(title, " – single-cell validation by sex"), outer=TRUE, font=2, cex=0.9)
}
for (co in COH) {
  x <- rd(file.path(SC, co$csv)); if (is.null(x)) next
  names(x)[names(x)=="n_Tcells"] <- "n_T"
  if (!is.null(co$sub) && names(co$sub) %in% names(x)) x <- x[as.character(x[[names(co$sub)]])==unname(co$sub),,drop=FALSE]
  if ("n_T" %in% names(x)) x <- x[!is.na(x$n_T) & x$n_T>=MIN_T,]
  x <- x[x$sex %in% c("M","F"),]
  cols <- names(RMAP)[names(RMAP) %in% names(x)]
  cols <- cols[vapply(cols, function(c) any(is.finite(x[[c]])), logical(1))]
  if (nrow(x) && length(cols)) save_fig(dotpanel(x, cols, co$id), paste0("SuppFig_S1_", co$id), 2.3*length(cols), 3.6)
}

## ---- S2: single-cell meta forest -------------------------------------------
eff <- rd(file.path(SC, "meta/sc_meta_effects.csv"))
if (!is.null(eff)) save_fig(function() {
  rds <- c("MAIT","Treg","CD8"); rds <- rds[rds %in% eff$readout]
  par(mfrow=c(1,length(rds)), mar=c(4,7,3,1), cex=0.8)
  for (rd_ in rds) {
    s <- eff[eff$readout==rd_ & is.finite(eff$delta),]; s <- s[order(s$cohort=="POOLED"),]
    n <- nrow(s); yy <- seq_len(n)
    plot(0,0,type="n", xlim=c(-1,1), ylim=c(0.5,n+0.5), yaxt="n", ylab="",
         xlab="Cliff's delta (>0 male)", main=rd_)
    abline(v=0, col="grey60")
    for (i in yy) { pooled <- s$cohort[i]=="POOLED"; col <- if(pooled)"#111111" else if(s$delta[i]>0)MALE else FEM
      segments(s$lo[i], i, s$hi[i], i, col=col, lwd=if(pooled)2.4 else 1.4)
      points(s$delta[i], i, pch=if(pooled)18 else 19, cex=if(pooled)1.6 else 1.1, col=col) }
    axis(2, at=yy, labels=s$cohort, las=1, cex.axis=0.7)
  }
}, "SuppFig_S2_sc_meta_forest", 11, 3.4)

## ---- S3: Guilliams QC (capture by sex + sex x diet) ------------------------
g <- rd(file.path(SC, "guilliams/out/Guilliams_clean_per_patient.csv"))
if (!is.null(g)) { g <- g[g$sex %in% c("M","F"),]
  save_fig(function() {
    par(mfrow=c(1,2), mar=c(4,4,3,1), cex=0.85)
    boxplot(n_T ~ sex, data=g, col=adjustcolor(c(FEM,MALE),.4), ylab="T cells captured per patient",
            main="A  Immune capture by sex", outline=FALSE)
    stripchart(n_T ~ sex, data=g, add=TRUE, vertical=TRUE, method="jitter", pch=19, cex=.7,
               col=c(FEM,MALE))
    if ("diet" %in% names(g)) { tb <- table(g$sex, g$diet)
      barplot(tb, beside=TRUE, col=c(FEM,MALE), legend.text=rownames(tb), xlab="diet",
              ylab="patients", main="B  Sex x diet (confound check)") }
    else plot.new()
  }, "SuppFig_S3_guilliams_QC", 8, 3.6)
}

cat("\nSupplementary figures written to", FIG, "as .pdf/.tiff/.png\n")
