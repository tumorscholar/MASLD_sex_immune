## Figure 4 — Disease-specificity vs disease-free GTEx liver + robustness (3.5)
## SELF-CONTAINED: Source in RStudio. Panels: A MASLD-vs-GTEx sex effect,
## B cross-cohort robustness (leave-one-cohort-out). The metabolic-deconfounding
## panel is now a separate SUPPLEMENTARY figure (suppfig_deconfounding.R).
## Reads gtex_healthy_sex.csv and maineffect_results.csv.
## Writes Figure_4_GTEx_specificity.pdf + .tiff (600 dpi).
## ---------------------------------------------------------------------------
RD  <- Sys.getenv("MASLD_REALDIR", "/data/Blizard-AlazawiLab/rk/MASLD_sex_meta")
FIG <- file.path(RD, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
MALE <- "#2C7FB8"; FEM <- "#D95F0E"; GREY <- "#888888"; GREEN <- "#4D9221"; AMB <- "#FDAE61"
pdf.options(family = "Helvetica")
save_fig <- function(draw, name, w, h) {
  pdf (file.path(FIG, paste0(name, ".pdf")),  width = w, height = h, family = "Helvetica"); draw(); dev.off()
  tiff(file.path(FIG, paste0(name, ".tiff")), width = w, height = h, units = "in", res = 600, compression = "lzw"); draw(); dev.off()
  cat("wrote", name, "-> pdf + tiff (600 dpi) in", FIG, "\n")
}
rd <- function(p) if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else { cat("MISSING:", p, "\n"); NULL }
LAB <- c(ct_MAIT = "MAIT", ct_MAITspec = "MAIT (SLC4A10/TRAV1-2)", ct_MAITpromisc = "MAIT (shared markers)",
  ct_Treg = "Treg", ct_CD8T = "CD8 T", ct_DC = "Dendritic", ct_NK = "NK", ct_Bcell = "B cell",
  ct_MonoMac = "Mono/Mac", st_Th1 = "Th1")
lab <- function(r) ifelse(r %in% names(LAB), LAB[r], r)
SHORT <- c(ct_MAIT = "MAIT", ct_MAITspec = "MAIT (recept.)", ct_MAITpromisc = "MAIT (shared)",
  ct_Treg = "Treg", ct_CD8T = "CD8 T", ct_DC = "cDC", ct_NK = "NK", ct_Bcell = "B", ct_MonoMac = "Mono/Mac", st_Th1 = "Th1")
slab <- function(r) ifelse(r %in% names(SHORT), SHORT[r], lab(r))
ROBUST <- "#5B7FA6"                                   # neutral (was green): fully robust
## dependency-free label de-overlap: iteratively push labels apart, draw leader lines
repel <- function(x, y, labels, cex = 0.6, col = "black", im = 300) {
  u <- par("usr"); dx <- u[2]-u[1]; dy <- u[4]-u[3]
  nx <- (x-u[1])/dx; ny <- (y-u[3])/dy; lx <- nx; ly <- ny
  wl <- strwidth(labels, cex=cex)/dx/2 + 0.006; hl <- strheight(labels, cex=cex)/dy/2 + 0.010
  for (it in 1:im) for (i in seq_along(lx)) {
    fx <- 0; fy <- 0
    for (j in seq_along(lx)) if (i!=j) {
      ax <- lx[i]-lx[j]; ay <- ly[i]-ly[j]
      ox <- (wl[i]+wl[j])-abs(ax); oy <- (hl[i]+hl[j])-abs(ay)
      if (ox>0 && oy>0) { fx <- fx + sign(ax+1e-6)*ox*0.45; fy <- fy + sign(ay+1e-6)*oy*0.45 }
    }
    lx[i] <- min(1-wl[i], max(wl[i], lx[i]+fx - (lx[i]-nx[i])*0.03))
    ly[i] <- min(1-hl[i], max(hl[i], ly[i]+fy - (ly[i]-ny[i])*0.03))
  }
  X <- u[1]+lx*dx; Y <- u[3]+ly*dy
  segments(x, y, X, Y, col="grey65", lwd=0.4); text(X, Y, labels, cex=cex, col=col)
}
HEAD <- c("ct_MAIT", "ct_MAITspec", "ct_Treg", "ct_CD8T", "ct_DC", "st_Th1")

me  <- rd(file.path(RD, "maineffect_results.csv")); stopifnot(!is.null(me))
meo <- me[me$coding == "fib_ord", ]; rownames(meo) <- meo$readout
gtx <- rd(file.path(RD, "gtex_healthy_sex.csv")); stopifnot(!is.null(gtx)); rownames(gtx) <- gtx$readout
NCOH <- max(meo$loco_n_sig, na.rm = TRUE); if (!is.finite(NCOH) || NCOH < 4) NCOH <- 6

save_fig(function() {
  layout(matrix(1:2, 1, 2), widths = c(1.25, 1)); par(mar = c(4, 4, 3, 1), cex = 0.72)
  ## A — MASLD vs GTEx
  common <- intersect(rownames(gtx), rownames(meo))
  xs <- meo[common, "beta_sexM"]; ys <- gtx[common, "beta"]; lim <- max(abs(c(xs, ys)), na.rm = TRUE) * 1.15
  plot(xs, ys, type = "n", xlim = c(-lim, lim), ylim = c(-lim, lim),
       xlab = "MASLD beta (>0 male)", ylab = "GTEx disease-free beta (>0 male)", main = "A  Disease-specificity")
  abline(0, 1, lty = 2, col = GREY); abline(h = 0, v = 0, col = "grey40", lwd = .5)
  hdi <- which(common %in% HEAD)
  for (i in seq_along(common)) { r <- common[i]; hd <- r %in% HEAD
    points(xs[i], ys[i], pch = 19, cex = if (hd) 1.4 else 0.6, col = if (hd) (if (xs[i] > 0) MALE else FEM) else adjustcolor(GREY, .5)) }
  repel(xs[hdi], ys[hdi], slab(common[hdi]), cex = 0.6)
  ## B — cross-cohort robustness (leave-one-cohort-out)
  hs <- HEAD[HEAD %in% rownames(meo)]; v <- meo[hs, "loco_n_sig"]
  bp <- barplot(rev(v), horiz = TRUE, col = rev(ifelse(v >= NCOH, ROBUST, AMB)), names.arg = rev(slab(hs)),
                las = 1, xlim = c(0, NCOH + 0.6), xlab = "cohorts significant (leave-one-out)",
                main = "B  Cross-cohort robustness", cex.names = 0.7)
  text(rev(v) + 0.1, bp, sprintf("%d/%d", rev(v), NCOH), adj = 0, cex = 0.75)
}, "Figure_4_GTEx_specificity", 7.1, 3.4)
