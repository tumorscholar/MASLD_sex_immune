## Figure 6 — Single-cell hepatic immune fractions by sex (Results 3.7) -------
## SELF-CONTAINED: Source in RStudio. Reads the in-house CITE-seq per-donor
## fractions CSV (set OWN_CITESEQ_CSV or MASLD_REALDIR). Writes
## Figure_6_singlecell.pdf + .tiff (600 dpi). One panel per readout.
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

## Input CSV: prefer an explicit OWN_CITESEQ_CSV, then the pipeline output (from
## sc_build_owncohort.R, RUN_SC), then the per-donor fractions file that SHIPS in
## the repo (data/own_cohort_percell_fractions.csv) so Fig 6 renders even when the
## single-cell steps were skipped.
cands <- c(Sys.getenv("OWN_CITESEQ_CSV", ""),
           file.path(RD, "single_cell/own_citeseq/own_per_donor_fractions.csv"),
           "../data/own_cohort_percell_fractions.csv",
           "data/own_cohort_percell_fractions.csv")
OWN_CSV <- cands[cands != "" & file.exists(cands)][1]
if (is.na(OWN_CSV)) stop("No in-house per-donor CSV found. Set OWN_CITESEQ_CSV, or ship data/own_cohort_percell_fractions.csv.")
cat("Fig 6 reading:", OWN_CSV, "\n")
own <- read.csv(OWN_CSV, stringsAsFactors = FALSE)

sexcol <- intersect(c("sex", "sex_assigned", "Sex"), names(own))[1]
## Print the actual in-house donor split for reference (Fig 6 / Results 3.7 = 5 M / 14 F).
if (!is.na(sexcol)) {
  .m <- sum(own[[sexcol]] %in% c("M","male","Male")); .f <- sum(own[[sexcol]] %in% c("F","female","Female"))
  cat(sprintf("[CHECK] in-house cohort donor split: %d male / %d female (n=%d)\n", .m, .f, nrow(own)))
}
## Accept either the pipeline names (MAIT_pctT ...) or the shipped provenance
## names (MAIT_ofT / cDC_ofImm ...). Values in 0-1 are displayed as %.
spec <- list(list(lab="MAIT (% T)",  alt=c("MAIT_pctT","MAIT_ofT")),
             list(lab="Treg (% T)",  alt=c("Treg_pctT","Treg_ofT")),
             list(lab="CD8 (% T)",   alt=c("CD8_pctT","CD8_ofT")),
             list(lab="cDC (% imm)", alt=c("cDC_pctImm","cDC_ofImm")))
resolved <- list()
for (s in spec) { col <- intersect(s$alt, names(own))[1]; if (!is.na(col)) resolved[[length(resolved)+1]] <- list(col=col, lab=s$lab) }
stopifnot(!is.na(sexcol), length(resolved) > 0)

getv <- function(col, idx) { v <- own[[col]][idx]; if (max(own[[col]], na.rm=TRUE) <= 1) v * 100 else v }  # 0-1 fraction -> %

save_fig(function() {
  par(mfrow = c(1, length(resolved)), mar = c(3, 4.3, 3.2, 0.6), mgp = c(2.6, 0.7, 0), oma = c(0, 0, 1.8, 0), cex = 0.72)
  isF <- own[[sexcol]] %in% c("F", "female", "Female"); isM <- own[[sexcol]] %in% c("M", "male", "Male")
  for (i in seq_along(resolved)) {
    rv <- resolved[[i]]$col; lab <- resolved[[i]]$lab
    f <- getv(rv, isF); m <- getv(rv, isM); p <- tryCatch(wilcox.test(m, f)$p.value, error = function(e) NA)
    plot(0, 0, type = "n", xlim = c(-0.6, 1.6), ylim = range(c(f, m), na.rm = TRUE), xaxt = "n", xlab = "",
         ylab = if (i == 1) "per-donor fraction (%)" else "", main = sprintf("%s\np=%.2g", lab, p))
    points(jitter(rep(0, length(f)), amount = .08), f, pch = 19, cex = .7, col = adjustcolor(FEM, .7))
    points(jitter(rep(1, length(m)), amount = .08), m, pch = 19, cex = .7, col = adjustcolor(MALE, .7))
    segments(-0.25, median(f, na.rm = TRUE), 0.25, median(f, na.rm = TRUE), lwd = 1.6)
    segments(0.75, median(m, na.rm = TRUE), 1.25, median(m, na.rm = TRUE), lwd = 1.6)
    axis(1, at = c(0, 1), labels = c("F", "M"))
  }
  mtext("In-house CITE-seq: hepatic immune fractions by sex", outer = TRUE, font = 2, cex = 0.85)
}, "Figure_6_singlecell", 2.0 * length(resolved), 3.2)
