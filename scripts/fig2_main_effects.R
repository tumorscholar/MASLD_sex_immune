## Figure 2 — Sex main effects: forest + per-sample distributions (Results 3.2)
## SELF-CONTAINED: Source in RStudio. Merges the old Fig 2 (forest) and old
## Fig 3 (distributions) into ONE multi-panel figure (A forest, B distributions).
## Reads maineffect_results.csv, maineffect_concordance.csv, analysis_matrix.csv.
## Writes Figure_2_main_effects.pdf (vector) + .tiff (600 dpi) to <REALDIR>/figures.
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
LAB <- c(ct_MAIT = "MAIT", ct_MAITspec = "MAIT (SLC4A10/TRAV1-2)", ct_MAITpromisc = "MAIT (shared markers)",
  ct_Treg = "Treg", ct_CD8T = "CD8 T", ct_DC = "Dendritic", ct_Tcell = "T cell", ct_NK = "NK", ct_Bcell = "B cell",
  ct_Plasma = "Plasma", ct_MonoMac = "Mono/Mac", ct_Neutrophil = "Neutrophil", st_Th1 = "Th1", st_Th17 = "Th17",
  st_CD8_cytotox = "CD8 cytotoxic", st_cytotoxCD4 = "Cytotoxic CD4", st_senescence = "Senescence",
  st_exhaustion = "Exhaustion", st_Tpex = "Tpex", st_Trm = "Trm")
lab <- function(r) ifelse(r %in% names(LAB), LAB[r], r)
SHORT <- c(ct_MAIT = "MAIT", ct_MAITspec = "MAIT (recept.)", ct_Treg = "Treg", ct_CD8T = "CD8 T", ct_DC = "cDC", st_Th1 = "Th1")
HEAD <- c("ct_MAIT", "ct_MAITspec", "ct_Treg", "ct_CD8T", "ct_DC", "st_Th1")

me  <- rd(file.path(RD, "maineffect_results.csv")); stopifnot(!is.null(me))
meo <- me[me$coding == "fib_ord", ]; rownames(meo) <- meo$readout
conc <- rd(file.path(RD, "maineffect_concordance.csv"))
mtx  <- rd(file.path(RD, "analysis_matrix.csv")); stopifnot(!is.null(mtx))
NSAMP <- nrow(mtx); NCOH <- length(unique(mtx$cohort))

save_fig(function() {
  ## row 1 = forest (spans 6 cols), row 2 = empty spacer, row 3 = 6 distribution panels
  layout(matrix(c(rep(1, 6), rep(0, 6), 2:7), 3, 6, byrow = TRUE), heights = c(2.2, 0.22, 1.05))
  ## ---- Panel A: forest across all readouts ----
  par(mar = c(5.6, 8.5, 3, 3), cex = 0.7)
  o <- meo[order(meo$beta_sexM), ]
  hd <- if (!is.null(conc)) conc$readout[conc$headline %in% c(TRUE, "True", "TRUE")] else HEAD
  n <- nrow(o); yy <- seq_len(n)
  xr <- range(c(o$beta_sexM - 1.96 * o$se, o$beta_sexM + 1.96 * o$se), na.rm = TRUE)
  plot(0, 0, type = "n", xlim = c(xr[1] - 0.34, xr[2] + 0.24), ylim = c(0.5, n + 0.5), yaxt = "n", ylab = "",
       xlab = "sex effect  beta  (< 0 female,  > 0 male)",
       main = sprintf("A  Pooled sex main effect across %d MASLD cohorts (n=%d, fibrosis-adjusted)", NCOH, NSAMP))
  abline(v = 0, lwd = 0.8)
  for (i in yy) { r <- rownames(o)[i]; b <- o$beta_sexM[i]; se <- o$se[i]; isH <- r %in% hd; col <- if (b > 0) MALE else FEM
    segments(b - 1.96 * se, i, b + 1.96 * se, i, col = col, lwd = if (isH) 2.4 else 1.1)
    points(b, i, pch = 19, cex = if (isH) 1.2 else 0.7, col = col)
    lc <- if (isH) sprintf("%s  LOCO %d/%d", star(o$fdr[i]), o$loco_n_sig[i], NCOH) else star(o$fdr[i])
    text(if (b > 0) b + 1.96 * se + 0.02 else b - 1.96 * se - 0.02, i, lc, adj = if (b > 0) 0 else 1, cex = 0.6) }
  axis(2, at = yy, labels = paste0(lab(rownames(o)), ifelse(rownames(o) %in% hd, " *", "")), las = 1, cex.axis = 0.72)
  mtext("* headline   *FDR<.05 **<.01 ***<.001   LOCO = leave-one-cohort-out fits significant", side = 1, line = 4.6, cex = 0.5, adj = 0)
  ## ---- Panel B: per-sample distributions of the headline readouts ----
  present <- HEAD[HEAD %in% names(mtx)]
  mz <- mtx; for (r in present) mz[[r]] <- ave(mz[[r]], mz$cohort, FUN = zc)
  par(mar = c(2.6, 3, 3.4, 0.4), cex = 0.58)
  for (k in seq_along(present)) {
    r <- present[k]
    f <- mz[[r]][mz$sex_assigned == "F"]; m <- mz[[r]][mz$sex_assigned == "M"]
    plot(0, 0, type = "n", xlim = c(-0.6, 1.6), ylim = range(c(f, m), na.rm = TRUE), xaxt = "n",
         xlab = "", ylab = if (k == 1) "within-cohort z" else "",
         main = sprintf("%s%s\nbeta=%+.2f %s", if (k == 1) "B  " else "", SHORT[r], meo[r, "beta_sexM"], star(meo[r, "fdr"])))
    abline(h = 0, lty = 3, col = "grey50")
    violin(f, 0, FEM); violin(m, 1, MALE)
    points(jitter(rep(0, length(f)), amount = .06), f, pch = 19, cex = .22, col = adjustcolor(FEM, .5))
    points(jitter(rep(1, length(m)), amount = .06), m, pch = 19, cex = .22, col = adjustcolor(MALE, .5))
    axis(1, at = c(0, 1), labels = c("F", "M"))
  }
}, "Figure_2_main_effects", 7.1, 8.6)
