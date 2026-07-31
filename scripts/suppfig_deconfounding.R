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
LAB <- c(ct_MAIT = "MAIT", ct_MAITspec = "MAIT (SLC4A10/TRAV1-2)", ct_MAITpromisc = "MAIT (shared markers)",
  ct_Treg = "Treg", ct_CD8T = "CD8 T", ct_DC = "Dendritic", st_Th1 = "Th1")
lab <- function(r) ifelse(r %in% names(LAB), LAB[r], r)

deconf <- rd(file.path(RD, "deconf_within_cohort.csv")); stopifnot(!is.null(deconf))

save_fig(function() {
  par(mar = c(4, 4, 3, 1), cex = 0.8)
  dh <- deconf[deconf$headline %in% c(TRUE, "True", "TRUE"), ]
  mx <- max(abs(c(dh$beta_unadj, dh$beta_adj)), na.rm = TRUE) * 1.2
  plot(dh$beta_unadj, dh$beta_adj, type = "n", xlim = c(-mx, mx), ylim = c(-mx, mx),
       xlab = "beta unadjusted", ylab = "beta + BMI/T2D/age", main = "Metabolic deconfounding")
  abline(0, 1, lty = 2, col = GREY); abline(h = 0, v = 0, col = "grey40", lwd = .4)
  points(dh$beta_unadj, dh$beta_adj, pch = 19, cex = 1.1, col = ifelse(dh$beta_unadj > 0, MALE, FEM))
  text(dh$beta_unadj, dh$beta_adj, lab(dh$readout), pos = 4, cex = 0.6, offset = 0.3)
}, "SuppFig_S4_metabolic_deconfounding", 3.6, 3.4)
