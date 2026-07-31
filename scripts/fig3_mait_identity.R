## Figure 3 — MAIT receptor-identity specificity (Results 3.3) ----------------
## SELF-CONTAINED: Source in RStudio. Reads analysis_matrix.csv and
## maineffect_results.csv. Writes Figure_3_MAIT_identity.pdf + .tiff (600 dpi).
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
star <- function(f) vapply(f, function(x) if (is.na(x)) "" else if (x < 1e-3) "***" else if (x < 1e-2) "**" else if (x < 0.05) "*" else if (x < 0.10) "ns" else "", character(1))
zc <- function(x) { m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) x * 0 else (x - m) / s }
violin <- function(x, at, col, width = 0.38) { x <- x[is.finite(x)]; if (length(x) < 3) return(invisible()); d <- density(x); xv <- d$y / max(d$y) * width; polygon(c(at - xv, rev(at + xv)), c(d$x, rev(d$x)), col = adjustcolor(col, 0.35), border = NA); segments(at - 0.28, median(x), at + 0.28, median(x), lwd = 1.4) }
SETLAB <- c(ct_MAITspec = "SLC4A10/TRAV1-2", ct_MAIT = "full panel", ct_MAITpromisc = "shared markers")

me  <- rd(file.path(RD, "maineffect_results.csv")); stopifnot(!is.null(me))
meo <- me[me$coding == "fib_ord", ]; rownames(meo) <- meo$readout
mtx <- rd(file.path(RD, "analysis_matrix.csv")); stopifnot(!is.null(mtx))

save_fig(function() {
  trio <- c("ct_MAITspec", "ct_MAIT", "ct_MAITpromisc"); trio <- trio[trio %in% names(mtx) & trio %in% rownames(meo)]
  layout(matrix(1:2, 1, 2), widths = c(1.3, 1)); par(mar = c(4, 4, 3.5, 1), cex = 0.72)
  mz <- mtx; for (t in trio) mz[[t]] <- ave(mz[[t]], mz$cohort, FUN = zc)
  pos <- 0; xt <- c(); xl <- c(); cen <- c()
  plot(0, 0, type = "n", xlim = c(-0.6, length(trio) * 2.7), ylim = range(unlist(mz[trio]), na.rm = TRUE) * 1.15,
       xaxt = "n", ylab = "within-cohort z-score", main = "A  MAIT signal is receptor-identity-specific")
  abline(h = 0, lty = 3, col = "grey50")
  for (t in trio) { cen <- c(cen, pos + 0.5)
    for (sc in c("F", "M")) { d <- mz[[t]][mz$sex_assigned == sc]; violin(d, pos, ifelse(sc == "M", MALE, FEM))
      xt <- c(xt, pos); xl <- c(xl, sc); pos <- pos + 1 }
    pos <- pos + 0.7 }
  axis(1, at = xt, labels = xl); top <- par("usr")[4]
  text(cen, top * 0.96, SETLAB[trio], font = 2, cex = 0.66)
  b <- meo[trio, "beta_sexM"]
  bp <- barplot(rev(b), horiz = TRUE, col = MALE, names.arg = rev(SETLAB[trio]), las = 1,
                xlim = c(0, max(b) * 1.45), xlab = "male-bias beta", main = "B  Effect size by marker set", cex.names = 0.7)
  text(rev(b) + 0.01, bp, sprintf("beta=%+.2f %s", rev(b), star(meo[trio, "fdr"])), adj = 0, cex = 0.7)
}, "Figure_3_MAIT_identity", 7.1, 3.3)
