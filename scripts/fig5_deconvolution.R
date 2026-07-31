## Figure 5 — Benchmarked deconvolution concordance (Results 3.6) -------------
## SELF-CONTAINED: Source in RStudio. Reads deconv_concordance.csv.
## Writes Figure_5_deconvolution.pdf + .tiff (600 dpi). Single-column width.
## ---------------------------------------------------------------------------
RD  <- Sys.getenv("MASLD_REALDIR", "/data/Blizard-AlazawiLab/rk/MASLD_sex_meta")
FIG <- file.path(RD, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
MALE <- "#2C7FB8"; FEM <- "#D95F0E"
pdf.options(family = "Helvetica")
save_fig <- function(draw, name, w, h) {
  pdf (file.path(FIG, paste0(name, ".pdf")),  width = w, height = h, family = "Helvetica"); draw(); dev.off()
  tiff(file.path(FIG, paste0(name, ".tiff")), width = w, height = h, units = "in", res = 600, compression = "lzw"); draw(); dev.off()
  cat("wrote", name, "-> pdf + tiff (600 dpi) in", FIG, "\n")
}
rd <- function(p) if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else { cat("MISSING:", p, "\n"); NULL }
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

dcon <- rd(file.path(RD, "deconv_concordance.csv")); stopifnot(!is.null(dcon))

save_fig(function() {
  x <- dcon[dcon$method == "xCell" & is.finite(dcon$singscore_beta) & is.finite(dcon$deconv_beta), ]
  par(mar = c(4, 4, 2.6, 1), cex = 0.78)
  lim <- max(0.6, max(abs(c(x$singscore_beta, x$deconv_beta)), na.rm = TRUE) * 1.2)
  plot(0, 0, type = "n", xlim = c(-lim, lim), ylim = c(-lim, lim), xlab = "singscore sex beta (>0 male)",
       ylab = "xCell sex beta (>0 male)", main = "Deconvolution reproduces singscore directions", cex.main = 0.92)
  rect(0, 0, lim, lim, col = "#EEF5EE", border = NA); rect(-lim, -lim, 0, 0, col = "#EEF5EE", border = NA)
  abline(h = 0, v = 0, lwd = .6)
  points(x$singscore_beta, x$deconv_beta, pch = 21, cex = 1.3, bg = ifelse(x$singscore_beta > 0, MALE, FEM))
  repel(x$singscore_beta, x$deconv_beta, x$deconv_celltype, cex = 0.58)
}, "Figure_5_deconvolution", 4.6, 4.2)
