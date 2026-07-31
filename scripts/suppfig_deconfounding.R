## Supplementary Figure S4 — Metabolic deconfounding (Results 3.4) ------------
## SELF-CONTAINED: Source in RStudio. This is Supplementary Figure S4 in the
## manuscript; §3.4 has no main figure (deconfounding is a robustness check).
## Reads deconf_within_cohort.csv.
## Writes SuppFig_S4_metabolic_deconfounding.pdf + .tiff (600 dpi).
## ---------------------------------------------------------------------------
RD  <- Sys.getenv("MASLD_REALDIR", "/data/Blizard-AlazawiLab/rk/MASLD_sex_meta")
FIG <- file.path(RD, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
MALE <- "#2C7FB8"; FEM <- "#D95F0E"; GREY <- "#888888"
pdf.options(family = "Helvetica")
save_fig <- function(draw, name, w, h) {
  pdf (file.path(FIG, paste0(name, ".pdf")),  width = w, height = h, family = "Helvetica"); draw(); dev.off()
  tiff(file.path(FIG, paste0(name, ".tiff")), width = w, height = h, units = "in", res = 600, compression = "lzw"); draw(); dev.off()
  cat("wrote", name, "-> pdf + tiff (600 dpi) in", FIG, "\n")
}
rd <- function(p) if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else { cat("MISSING:", p, "\n"); NULL }
LAB <- c(ct_MAIT = "MAIT", ct_MAITspec = "MAIT (recept.)", ct_MAITpromisc = "MAIT (shared)",
  ct_Treg = "Treg", ct_CD8T = "CD8 T", ct_DC = "Dendritic", ct_MonoMac = "Mono/Mac",
  st_Th1 = "Th1", st_exhaustion = "Exhaustion", st_Tpex = "Tpex")
lab <- function(r) ifelse(r %in% names(LAB), as.character(LAB[r]), r)
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

deconf <- rd(file.path(RD, "deconf_within_cohort.csv")); stopifnot(!is.null(deconf))

save_fig(function() {
  par(mar = c(4, 4, 3, 1), cex = 0.8)
  dh <- deconf[deconf$headline %in% c(TRUE, "True", "TRUE"), ]
  mx <- max(abs(c(dh$beta_unadj, dh$beta_adj)), na.rm = TRUE) * 1.2
  plot(dh$beta_unadj, dh$beta_adj, type = "n", xlim = c(-mx, mx), ylim = c(-mx, mx),
       xlab = "beta unadjusted", ylab = "beta + BMI/T2D/age", main = "Metabolic deconfounding")
  abline(0, 1, lty = 2, col = GREY); abline(h = 0, v = 0, col = "grey40", lwd = .4)
  points(dh$beta_unadj, dh$beta_adj, pch = 19, cex = 1.1, col = ifelse(dh$beta_unadj > 0, MALE, FEM))
  repel(dh$beta_unadj, dh$beta_adj, lab(dh$readout), cex = 0.6)
}, "SuppFig_S4_metabolic_deconfounding", 4.0, 3.8)
