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

dcon <- rd(file.path(RD, "deconv_concordance.csv")); stopifnot(!is.null(dcon))

save_fig(function() {
  x <- dcon[dcon$method == "xCell" & is.finite(dcon$singscore_beta) & is.finite(dcon$deconv_beta), ]
  par(mar = c(4, 4, 3, 1), cex = 0.8)
  lim <- max(0.6, max(abs(c(x$singscore_beta, x$deconv_beta)), na.rm = TRUE) * 1.2)
  plot(0, 0, type = "n", xlim = c(-lim, lim), ylim = c(-lim, lim), xlab = "singscore sex beta (>0 male)",
       ylab = "xCell sex beta (>0 male)", main = "Deconvolution reproduces singscore sex directions")
  rect(0, 0, lim, lim, col = "#EEF5EE", border = NA); rect(-lim, -lim, 0, 0, col = "#EEF5EE", border = NA)
  abline(h = 0, v = 0, lwd = .6)
  points(x$singscore_beta, x$deconv_beta, pch = 21, cex = 1.4, bg = ifelse(x$singscore_beta > 0, MALE, FEM))
  text(x$singscore_beta, x$deconv_beta, x$deconv_celltype, pos = 4, cex = 0.6, offset = 0.3)
}, "Figure_5_deconvolution", 3.5, 3.3)
