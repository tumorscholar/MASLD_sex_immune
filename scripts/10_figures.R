## 10_figures.R ---------------------------------------------------------------
## Regenerates the main-text figures from the pipeline output CSVs and saves each
## as PDF (vector) + TIFF (300 dpi, LZW) + PNG – the formats needed for journal
## submission. Run after the analysis pipeline (03/04/08) has written its CSVs.
##   Fig 1  expression-based sex assignment      (per_sample_sex_calls.csv)
##   Fig 2  pooled sex main-effect forest        (maineffect_results/-concordance)
##   Fig 3  per-sample distributions             (analysis_matrix.csv)
##   Fig 4  MAIT receptor-identity specificity   (analysis_matrix + maineffect)
##   Fig 5  disease-specificity + robustness     (gtex_healthy_sex, deconf, maineffect)
##   Fig 6  deconvolution concordance            (deconv_concordance.csv)
## Main Figure 7 (combined in-house + public single-cell dot plot) is produced by
## sc_meta_dotplot.R. This script also regenerates the in-house panel of
## Supplementary Figure S1 (see the S1 block below).
## ===========================================================================
RD     <- Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta")
SEXDIR <- if (dir.exists(file.path(RD,"sex_assign_out"))) file.path(RD,"sex_assign_out") else path.expand("~/sex_assign_out")
FIG    <- file.path(RD, "figures"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
MALE <- "#2c7fb8"; FEM <- "#d95f0e"; GREY <- "#888888"

## ---- save each figure in three formats -------------------------------------
save_fig <- function(draw, name, w, h) {
  render <- function() if (inherits(draw,"ggplot")) print(draw) else draw()
  pdf (file.path(FIG, paste0(name,".pdf")),  width=w, height=h);                                   render(); dev.off()
  tiff(file.path(FIG, paste0(name,".tiff")), width=w, height=h, units="in", res=300, compression="lzw"); render(); dev.off()
  png (file.path(FIG, paste0(name,".png")),  width=w, height=h, units="in", res=300);              render(); dev.off()
  cat("wrote", name, "(pdf, tiff, png)\n")
}
star <- function(f) vapply(f, function(x) if(is.na(x)) "" else if(x<1e-3)"***" else if(x<1e-2)"**" else if(x<0.05)"*" else if(x<0.10)"ns" else "", character(1))
LAB <- c(ct_MAIT="MAIT", ct_MAITspec="MAIT (SLC4A10/TRAV1-2)", ct_MAITpromisc="MAIT (promiscuous)",
  ct_Treg="Treg", ct_CD8T="CD8 T", ct_DC="Dendritic", ct_Tcell="T cell", ct_NK="NK", ct_Bcell="B cell",
  ct_Plasma="Plasma", ct_MonoMac="Mono/Mac", ct_Neutrophil="Neutrophil", st_Th1="Th1", st_Th17="Th17",
  st_CD8_cytotox="CD8 cytotoxic", st_cytotoxCD4="Cytotoxic CD4", st_senescence="Senescence",
  st_exhaustion="Exhaustion", st_Tpex="Tpex", st_Trm="Trm")
lab <- function(r) ifelse(r %in% names(LAB), LAB[r], r)
HEAD <- c("ct_MAIT","ct_MAITspec","ct_Treg","ct_CD8T","ct_DC","st_Th1")
zc <- function(x) { m<-mean(x,na.rm=TRUE); s<-sd(x,na.rm=TRUE); if(is.na(s)||s==0) x*0 else (x-m)/s }
violin <- function(x, at, col, width=0.38) {
  x <- x[is.finite(x)]; if (length(x) < 3) return(invisible())
  d <- density(x); xv <- d$y/max(d$y)*width
  polygon(c(at-xv, rev(at+xv)), c(d$x, rev(d$x)), col=adjustcolor(col,0.35), border=NA)
  segments(at-0.28, median(x), at+0.28, median(x), lwd=1.4)
}

## ---- read result CSVs ------------------------------------------------------
rd <- function(p) if (file.exists(p)) read.csv(p, stringsAsFactors=FALSE) else { cat("missing:",p,"\n"); NULL }
calls  <- rd(file.path(SEXDIR,"per_sample_sex_calls.csv"))
mtx    <- rd(file.path(RD,"analysis_matrix.csv"))
me     <- rd(file.path(RD,"maineffect_results.csv"))
meo    <- if(!is.null(me)) { m<-me[me$coding=="fib_ord",]; rownames(m)<-m$readout; m } else NULL
conc   <- rd(file.path(RD,"maineffect_concordance.csv"))
deconf <- rd(file.path(RD,"deconf_within_cohort.csv"))
gtx    <- rd(file.path(RD,"gtex_healthy_sex.csv")); if(!is.null(gtx)) rownames(gtx)<-gtx$readout
dcon   <- rd(file.path(RD,"deconv_concordance.csv"))

## ===================== FIG 1  sex assignment =====================
if (!is.null(calls)) save_fig(function() {
  layout(matrix(1:3,1,3), widths=c(1.25,1,1)); par(mar=c(4,4,3,1), cex=0.85)
  ## A scatter
  plot(calls$zXIST, calls$zYmean, type="n", xlab="XIST z-score (within cohort)",
       ylab="Y-panel mean z-score", main="A  Expression-based sex calls")
  abline(h=0,v=0,lty=3,col="grey40")
  for (sc in c("M","F")) { d<-calls[calls$sex_assigned==sc,]
    points(d$zXIST,d$zYmean,pch=19,cex=0.5,col=adjustcolor(ifelse(sc=="M",MALE,FEM),0.7)) }
  legend("topright", pch=19, col=c(MALE,FEM), bty="n", cex=0.8,
         legend=c(paste0("M (n=",sum(calls$sex_assigned=="M"),")"), paste0("F (n=",sum(calls$sex_assigned=="F"),")")))
  ## B concordance per cohort
  cc <- calls[calls$recorded_sex %in% c("M","F"),]
  if (nrow(cc)) {
    g <- tapply(cc$sex_assigned==cc$recorded_sex, cc$gse, function(x) 100*mean(x))
    n <- tapply(cc$gse, cc$gse, length)
    bp <- barplot(g, horiz=TRUE, col="#4d9221", xlim=c(0,108), names.arg=paste0(names(g),"\n(n=",n[names(g)],")"),
                  las=1, cex.names=0.6, xlab="concordance with recorded sex (%)", main="B  Validation vs metadata")
    text(pmin(g,99)+2, bp, sprintf("%.0f%%",g), cex=0.7, adj=0)
  } else { plot.new(); title("B  (no recorded sex)") }
  ## C separation histogram
  sep <- calls$zYmean - calls$zXIST
  hM <- hist(sep[calls$sex_assigned=="M"], breaks=25, plot=FALSE)
  hF <- hist(sep[calls$sex_assigned=="F"], breaks=25, plot=FALSE)
  plot(0,0,type="n", xlim=range(sep,na.rm=TRUE), ylim=c(0,max(hM$counts,hF$counts)),
       xlab="Y-panel minus XIST (z)", ylab="samples", main="C  Bimodal separation")
  plot(hF, add=TRUE, col=adjustcolor(FEM,0.6), border=NA); plot(hM, add=TRUE, col=adjustcolor(MALE,0.6), border=NA)
  abline(v=0,lty=2); legend("topright", fill=c(MALE,FEM), legend=c("M","F"), bty="n", cex=0.8)
}, "Fig1_sex_assignment", 11, 3.6)

## ===================== FIG 2  main-effect forest =====================
if (!is.null(meo)) save_fig(function() {
  o <- meo[order(meo$beta_sexM),]; hd <- if(!is.null(conc)) conc$readout[conc$headline %in% c(TRUE,"True","TRUE")] else HEAD
  par(mar=c(4,9,3,3), cex=0.85); n<-nrow(o); yy<-seq_len(n)
  xr <- range(c(o$beta_sexM-1.96*o$se, o$beta_sexM+1.96*o$se), na.rm=TRUE)
  plot(0,0,type="n", xlim=c(xr[1]-0.34,xr[2]+0.22), ylim=c(0.5,n+0.5), yaxt="n", ylab="",
       xlab="sex effect  beta  (<0 higher in FEMALES     >0 higher in MALES)",
       main="Pooled sex main effect across 6 MASLD cohorts (n=616, fibrosis-adjusted)")
  abline(v=0, lwd=0.8)
  for (i in yy) { r<-rownames(o)[i]; b<-o$beta_sexM[i]; se<-o$se[i]; isH<-r%in%hd; col<-if(b>0)MALE else FEM
    segments(b-1.96*se, i, b+1.96*se, i, col=col, lwd=if(isH)2.4 else 1.1)
    points(b, i, pch=19, cex=if(isH)1.3 else 0.7, col=col)
    lc <- if(isH) sprintf("%s  LOCO %d/4", star(o$fdr[i]), o$loco_n_sig[i]) else star(o$fdr[i])
    text(if(b>0) b+1.96*se+0.02 else b-1.96*se-0.02, i, lc, adj=if(b>0)0 else 1, cex=0.62) }
  axis(2, at=yy, labels=paste0(lab(rownames(o)), ifelse(rownames(o)%in%hd," *","")), las=1, cex.axis=0.7)
  mtext("* headline   *FDR<.05 **<.01 ***<.001   LOCO = leave-one-cohort-out", side=1, line=2.6, cex=0.55, adj=0)
}, "Fig2_maineffect_forest", 7.4, 6.2)

## ===================== FIG 3  distributions =====================
if (!is.null(mtx) && !is.null(meo)) save_fig(function() {
  present <- HEAD[HEAD %in% names(mtx)]
  mz <- mtx; for (r in present) mz[[r]] <- ave(mz[[r]], mz$cohort, FUN=zc)
  par(mfrow=c(1,length(present)), mar=c(3,3.2,3.5,0.6), oma=c(0,0,2,0), cex=0.8)
  for (r in present) {
    f <- mz[[r]][mz$sex_assigned=="F"]; m <- mz[[r]][mz$sex_assigned=="M"]
    yr <- range(c(f,m), na.rm=TRUE)
    plot(0,0,type="n", xlim=c(-0.6,1.6), ylim=yr, xaxt="n", xlab="", ylab=if(r==present[1])"within-cohort z-score" else "",
         main=sprintf("%s\nbeta=%+.2f %s", lab(r), meo[r,"beta_sexM"], star(meo[r,"fdr"])))
    abline(h=0, lty=3, col="grey50")
    violin(f,0,FEM); violin(m,1,MALE)
    points(jitter(rep(0,length(f)),amount=.06), f, pch=19, cex=.25, col=adjustcolor(FEM,.5))
    points(jitter(rep(1,length(m)),amount=.06), m, pch=19, cex=.25, col=adjustcolor(MALE,.5))
    axis(1, at=c(0,1), labels=c("F","M"))
  }
  mtext("Per-sample distribution of headline readouts by assigned sex", outer=TRUE, font=2, cex=0.9)
}, "Fig3_distributions", 12, 3.6)

## ===================== FIG 4  MAIT identity =====================
if (!is.null(mtx) && !is.null(meo)) save_fig(function() {
  trio <- c("ct_MAITspec","ct_MAIT","ct_MAITpromisc"); trio <- trio[trio %in% names(mtx) & trio %in% rownames(meo)]
  layout(matrix(1:2,1,2), widths=c(1.3,1)); par(mar=c(4,4,3.5,1), cex=0.85)
  mz <- mtx; for (t in trio) mz[[t]] <- ave(mz[[t]], mz$cohort, FUN=zc)
  pos<-0; xt<-c(); xl<-c(); cen<-c()
  plot(0,0,type="n", xlim=c(-0.6, length(trio)*2.7), ylim=range(unlist(mz[trio]),na.rm=TRUE)*1.15,
       xaxt="n", ylab="within-cohort z-score", main="A  MAIT signal is receptor-identity-specific")
  abline(h=0, lty=3, col="grey50")
  for (t in trio) { cen<-c(cen,pos+0.5)
    for (sc in c("F","M")) { d<-mz[[t]][mz$sex_assigned==sc]; violin(d,pos,ifelse(sc=="M",MALE,FEM))
      xt<-c(xt,pos); xl<-c(xl,sc); pos<-pos+1 }
    pos<-pos+0.7 }
  axis(1, at=xt, labels=xl); top<-par("usr")[4]
  text(cen, top*0.96, gsub("MAIT ","",lab(trio)), font=2, cex=0.7)
  b <- meo[trio,"beta_sexM"]
  bp <- barplot(rev(b), horiz=TRUE, col=MALE, names.arg=rev(gsub("MAIT ","",lab(trio))), las=1,
                xlim=c(0,max(b)*1.4), xlab="male-bias beta", main="B  Effect size by marker set", cex.names=0.7)
  text(rev(b)+0.01, bp, sprintf("beta=%+.2f %s", rev(b), star(meo[trio,"fdr"])), adj=0, cex=0.7)
}, "Fig4_MAIT_identity", 9, 3.7)

## ===================== FIG 5  disease-specificity + robustness =====================
if (!is.null(gtx) && !is.null(meo)) save_fig(function() {
  layout(matrix(1:3,1,3), widths=c(1.25,1,1)); par(mar=c(4,4,3,1), cex=0.8)
  ## A MASLD vs GTEx
  common <- intersect(rownames(gtx), rownames(meo))
  xs <- meo[common,"beta_sexM"]; ys <- gtx[common,"beta"]; lim<-max(abs(c(xs,ys)),na.rm=TRUE)*1.15
  plot(xs, ys, type="n", xlim=c(-lim,lim), ylim=c(-lim,lim),
       xlab="MASLD beta (>0 male)", ylab="GTEx healthy beta (>0 male)", main="A  Disease-specificity")
  abline(0,1,lty=2,col=GREY); abline(h=0,v=0,col="grey40",lwd=.5)
  for (i in seq_along(common)) { r<-common[i]; hd<-r%in%HEAD
    points(xs[i],ys[i], pch=19, cex=if(hd)1.4 else 0.6, col=if(hd)(if(xs[i]>0)MALE else FEM) else adjustcolor(GREY,.5))
    if (hd) text(xs[i],ys[i], lab(r), pos=4, cex=0.6, offset=0.3) }
  ## B deconfound
  if (!is.null(deconf)) { dh<-deconf[deconf$headline %in% c(TRUE,"True","TRUE"),]
    mx<-max(abs(c(dh$beta_unadj,dh$beta_adj)),na.rm=TRUE)*1.2
    plot(dh$beta_unadj, dh$beta_adj, type="n", xlim=c(-mx,mx), ylim=c(-mx,mx),
         xlab="beta unadjusted", ylab="beta + BMI/T2D/age", main="B  Metabolic deconfounding")
    abline(0,1,lty=2,col=GREY); abline(h=0,v=0,col="grey40",lwd=.4)
    points(dh$beta_unadj, dh$beta_adj, pch=19, cex=1.1, col=ifelse(dh$beta_unadj>0,MALE,FEM))
    text(dh$beta_unadj, dh$beta_adj, lab(dh$readout), pos=4, cex=0.55, offset=0.3)
  } else plot.new()
  ## C LOCO
  hs <- HEAD[HEAD %in% rownames(meo)]; v <- meo[hs,"loco_n_sig"]
  bp <- barplot(rev(v), horiz=TRUE, col=rev(ifelse(v>=3,"#4d9221","#fdae61")), names.arg=rev(lab(hs)),
                las=1, xlim=c(0,4.6), xlab="cohorts significant (leave-one-out)", main="C  Cross-cohort robustness", cex.names=0.7)
  text(rev(v)+0.1, bp, sprintf("%d/4",rev(v)), adj=0, cex=0.75)
}, "Fig5_disease_specificity", 12, 3.9)

## ===================== FIG 6  deconvolution concordance =====================
if (!is.null(dcon)) save_fig(function() {
  x <- dcon[dcon$method=="xCell" & is.finite(dcon$singscore_beta) & is.finite(dcon$deconv_beta), ]
  par(mar=c(4,4,3,1), cex=0.85)
  lim <- max(0.6, max(abs(c(x$singscore_beta,x$deconv_beta)),na.rm=TRUE)*1.2)
  plot(0,0,type="n", xlim=c(-lim,lim), ylim=c(-lim,lim), xlab="singscore sex beta (>0 male)",
       ylab="xCell sex beta (>0 male)", main="Benchmarked deconvolution reproduces singscore sex directions")
  rect(0,0,lim,lim,col="#eef5ee",border=NA); rect(-lim,-lim,0,0,col="#eef5ee",border=NA)
  abline(h=0,v=0,lwd=.6)
  points(x$singscore_beta, x$deconv_beta, pch=21, cex=1.4, bg=ifelse(x$singscore_beta>0,MALE,FEM))
  text(x$singscore_beta, x$deconv_beta, x$deconv_celltype, pos=4, cex=0.6, offset=0.3)
}, "Fig6_deconv_concordance", 5.6, 5.2)

## ============== SuppFig S1 (in-house panel)  CITE-seq validation ==============
## The MAIN Figure 7 (in-house + four public cohorts, one combined dot plot) is
## produced by sc_meta_dotplot.R. This block regenerates the in-house cohort on
## its own, as the in-house panel of Supplementary Figure S1 (per-cohort breakout;
## the four public panels come from 11_supp_figures.R). Reads the per-patient
## fraction table written by sc_build_owncohort.R (proportions shown as %).
## Default looks first for the shipped data/ copy, then the pipeline output.
OWN_CANDIDATES <- c(
  Sys.getenv("OWN_CITESEQ_CSV", ""),
  file.path(RD, "single_cell/own_citeseq/own_cohort_percell_fractions.csv"),
  file.path(dirname(RD), "data/own_cohort_percell_fractions.csv"),
  "../data/own_cohort_percell_fractions.csv",
  "data/own_cohort_percell_fractions.csv")
OWN_CSV <- OWN_CANDIDATES[OWN_CANDIDATES != "" & file.exists(OWN_CANDIDATES)][1]
own <- if (length(OWN_CSV) && !is.na(OWN_CSV)) rd(OWN_CSV) else NULL

if (!is.null(own)) {
  suppressMessages(library(ggplot2))
  sexcol <- intersect(c("sex","sex_assigned","Sex"), names(own))[1]
  ## accept either proportion (_ofT/_ofImm) or percent (_pctT/_pctImm) columns
  READ <- list(
    MAIT = c("MAIT (% T)",   c("MAIT_ofT","MAIT_pctT")),
    Treg = c("Treg (% T)",   c("Treg_ofT","Treg_pctT")),
    CD8  = c("CD8 (% T)",    c("CD8_ofT","CD8_pctT")),
    cDC  = c("cDC (% imm)",  c("cDC_ofImm","cDC_pctImm")))
  pick <- function(alts) { h <- alts[alts %in% names(own)][1]; if (is.na(h)) NULL else h }
  long <- list(); stats <- list()
  for (nm in names(READ)) {
    col <- pick(READ[[nm]][-1]); if (is.null(col)) next
    lab <- READ[[nm]][1]
    v <- own[[col]]; if (max(v, na.rm=TRUE) <= 1.5) v <- v * 100   # proportion -> %
    isF <- own[[sexcol]] %in% c("F","female","Female"); isM <- own[[sexcol]] %in% c("M","male","Male")
    long[[nm]] <- data.frame(readout=lab, sex=ifelse(isM,"M",ifelse(isF,"F",NA)), value=v)
    p <- tryCatch(wilcox.test(v[isM], v[isF])$p.value, error=function(e) NA)
    stats[[nm]] <- data.frame(readout=lab, p=p)
  }
  D <- do.call(rbind, long); D <- D[D$sex %in% c("M","F") & is.finite(D$value), ]
  S <- do.call(rbind, stats)
  lev <- vapply(READ, function(z) z[1], character(1)); lev <- lev[lev %in% D$readout]
  D$readout <- factor(D$readout, levels=lev); S$readout <- factor(S$readout, levels=lev)
  S$lab <- ifelse(is.na(S$p), "", paste0("p=", signif(S$p,2)))
  ytop <- aggregate(value~readout, D, function(v) max(v)*1.05); S <- merge(S, ytop, all.x=TRUE)
  nM <- sum(own[[sexcol]] %in% c("M","male","Male")); nF <- sum(own[[sexcol]] %in% c("F","female","Female"))

  p7 <- ggplot(D, aes(sex, value, colour=sex)) +
    geom_jitter(width=0.18, height=0, size=1.8, alpha=0.8) +
    stat_summary(fun=median, geom="crossbar", width=0.55, linewidth=0.3, colour="black") +
    geom_text(data=S, aes(x=1.5, y=value, label=lab), inherit.aes=FALSE, size=2.7, vjust=0, colour="grey25") +
    facet_wrap(~readout, scales="free_y", nrow=1) +
    scale_colour_manual(values=c(M=MALE, F=FEM)) +
    labs(x=NULL, y="per-patient fraction (%)",
         title="In-house CITE-seq cohort (Supplementary Figure S1 panel)",
         subtitle=sprintf("each dot = one patient (%d M / %d F) · black bar = median · blue male / orange female", nM, nF)) +
    theme_bw(base_size=9) +
    theme(legend.position="none", panel.grid.minor=element_blank(),
          strip.background=element_rect(fill="grey95", colour=NA),
          plot.subtitle=element_text(size=8, colour="grey40"))
  save_fig(p7, "SuppFig_S1_InHouse", 9, 3.4)
} else cat("In-house S1 panel skipped: no in-house fractions CSV found. Run sc_build_owncohort.R",
           "or set OWN_CITESEQ_CSV.\n")

cat("\nAll available figures written to", FIG, "as .pdf/.tiff/.png\n")
