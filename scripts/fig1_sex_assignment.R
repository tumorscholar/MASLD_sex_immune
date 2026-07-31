## Figure 1 — Expression-based sex assignment (Results 3.1) -------------------
## SELF-CONTAINED: open in RStudio and Source (or Run All). No other script
## needs to be run first. Reads per_sample_sex_calls.csv and writes
## Figure_1_sex_assignment.pdf (vector) + .tiff (600 dpi, LZW) to
## <MASLD_REALDIR>/figures. JHEP style: RGB, Helvetica/Arial, colourblind-safe
## blue = male, orange = female.
## ---------------------------------------------------------------------------
RD     <- Sys.getenv("MASLD_REALDIR", "/data/Blizard-AlazawiLab/rk/MASLD_sex_meta")
SEXDIR <- if (dir.exists(file.path(RD, "sex_assign_out"))) file.path(RD, "sex_assign_out") else path.expand("~/sex_assign_out")
FIG    <- file.path(RD, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
MALE <- "#2C7FB8"; FEM <- "#D95F0E"; GREEN <- "#4D9221"
pdf.options(family = "Helvetica")
save_fig <- function(draw, name, w, h) {
  pdf (file.path(FIG, paste0(name, ".pdf")),  width = w, height = h, family = "Helvetica"); draw(); dev.off()
  tiff(file.path(FIG, paste0(name, ".tiff")), width = w, height = h, units = "in", res = 600, compression = "lzw"); draw(); dev.off()
  cat("wrote", name, "-> pdf + tiff (600 dpi) in", FIG, "\n")
}
rd <- function(p) if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else { cat("MISSING:", p, "\n"); NULL }

calls <- rd(file.path(SEXDIR, "per_sample_sex_calls.csv"))
stopifnot(!is.null(calls))

save_fig(function() {
  layout(matrix(1:3, 1, 3), widths = c(1.25, 1, 1)); par(mar = c(4, 4, 3, 1), cex = 0.85)
  ## A — expression-based sex calls
  plot(calls$zXIST, calls$zYmean, type = "n", xlab = "XIST z-score (within cohort)",
       ylab = "Y-panel mean z-score", main = "A  Expression-based sex calls")
  abline(h = 0, v = 0, lty = 3, col = "grey40")
  for (sc in c("M", "F")) { d <- calls[calls$sex_assigned == sc, ]
    points(d$zXIST, d$zYmean, pch = 19, cex = 0.5, col = adjustcolor(ifelse(sc == "M", MALE, FEM), 0.7)) }
  legend("topright", pch = 19, col = c(MALE, FEM), bty = "n", cex = 0.8,
         legend = c(paste0("M (n=", sum(calls$sex_assigned == "M"), ")"), paste0("F (n=", sum(calls$sex_assigned == "F"), ")")))
  ## B — concordance with recorded sex, per cohort
  cc <- calls[calls$recorded_sex %in% c("M", "F"), ]
  if (nrow(cc)) {
    g <- tapply(cc$sex_assigned == cc$recorded_sex, cc$gse, function(x) 100 * mean(x))
    n <- tapply(cc$gse, cc$gse, length)
    bp <- barplot(g, horiz = TRUE, col = GREEN, xlim = c(0, 108), names.arg = paste0(names(g), "\n(n=", n[names(g)], ")"),
                  las = 1, cex.names = 0.6, xlab = "concordance with recorded sex (%)", main = "B  Validation vs metadata")
    text(pmin(g, 99) + 2, bp, sprintf("%.0f%%", g), cex = 0.7, adj = 0)
  } else { plot.new(); title("B  (no recorded sex)") }
  ## C — bimodal separation
  sep <- calls$zYmean - calls$zXIST
  hM <- hist(sep[calls$sex_assigned == "M"], breaks = 25, plot = FALSE)
  hF <- hist(sep[calls$sex_assigned == "F"], breaks = 25, plot = FALSE)
  plot(0, 0, type = "n", xlim = range(sep, na.rm = TRUE), ylim = c(0, max(hM$counts, hF$counts)),
       xlab = "Y-panel minus XIST (z)", ylab = "samples", main = "C  Bimodal separation")
  plot(hF, add = TRUE, col = adjustcolor(FEM, 0.6), border = NA); plot(hM, add = TRUE, col = adjustcolor(MALE, 0.6), border = NA)
  abline(v = 0, lty = 2); legend("topright", fill = c(MALE, FEM), legend = c("M", "F"), bty = "n", cex = 0.8)
}, "Figure_1_sex_assignment", 11, 3.6)
