## make_workflow_schematic.R --------------------------------------------------
## Draws Supplementary Figure S7 (analysis workflow) as PDF + TIFF + PNG using
## base graphics only (no extra packages). This is a schematic, not a data
## figure, so it takes no pipeline input. Layout mirrors the two-axis design:
## public data + GTEx -> sex assignment -> scoring -> mixed model -> meta ->
## four specificity/validation layers -> GTEx classification -> the two axes.
## ---------------------------------------------------------------------------
source("00_config.R")   # for REALDIR / save path convention (optional)

outdir <- tryCatch(file.path(REALDIR, "figures"), error = function(e) ".")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## palette (matches the manuscript table shades)
BLUE   <- "#DCE9F5"; BLUE_E   <- "#5A7DA0"
ORANGE <- "#FBE6D5"; ORANGE_E <- "#C08040"
GREEN  <- "#DDEBDD"; GREEN_E  <- "#5A8A5A"
GREY   <- "#F1F1F1"; GREY_E   <- "#8A8A8A"
NEUT   <- "#EAF0F6"; NEUT_E   <- "#6B7B8C"

draw <- function() {
  op <- par(mar = c(0,0,0,0)); on.exit(par(op))
  plot(NA, xlim = c(0,100), ylim = c(0,100), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")

  box <- function(xc, yc, w, h, fill, edge, head, sub = NULL, cexh = 0.95, cexs = 0.82) {
    rect(xc-w/2, yc-h/2, xc+w/2, yc+h/2, col = fill, border = edge, lwd = 1.6)
    if (is.null(sub)) {
      text(xc, yc, head, font = 2, cex = cexh)
    } else {
      text(xc, yc + h*0.16, head, font = 2, cex = cexh)
      text(xc, yc - h*0.20, sub, cex = cexs, col = "#333333")
    }
  }
  arr <- function(x1,y1,x2,y2, col = "#444444") arrows(x1,y1,x2,y2, length = 0.08, lwd = 1.8, col = col)

  ## inputs
  box(32, 93, 40, 8, BLUE,  BLUE_E,  "6 public MASLD liver datasets", "637 samples (616 QC); RNA-seq + array")
  box(76, 93, 34, 8, GREEN, GREEN_E, "GTEx v8 disease-free liver",    "226 donors (161 M / 65 F)")
  ## pipeline spine
  box(50, 80, 46, 7, GREY, GREY_E, "Expression-based sex assignment", "XIST vs Y-chromosome gene panel")
  arr(32,89,45,84); arr(76,89,60,84)
  box(50, 69, 46, 7, GREY, GREY_E, "Per-sample signature scoring (singscore)", "10 immune cell types + 8 states")
  arr(50,76.5,50,72.5)
  box(50, 58, 46, 7, GREY, GREY_E, "Sex effect per dataset", "fibrosis-adjusted mixed models")
  arr(50,65.5,50,61.5)
  box(50, 47, 46, 7.5, NEUT, NEUT_E, "Random-effects meta-analysis", "pooled beta, 95% CI, I2, prediction interval")
  arr(50,54.5,50,51)

  ## validation band
  vy <- 34; vw <- 21; xs <- c(14, 38, 62, 86)
  labs <- list(c("Parent-lineage","adjustment"), c("Metabolic","deconfounding"),
               c("Benchmarked","deconvolution"), c("Single-cell cohorts","in-house + 4 public"))
  arr(50,43,50,40.5)
  segments(xs[1], 40.5, xs[4], 40.5, lwd = 1.8, col = "#444444")
  for (i in seq_along(xs)) {
    segments(xs[i],40.5,xs[i],38.2, lwd = 1.8, col = "#444444"); arrows(xs[i],38.6,xs[i],vy+4, length=0.07,lwd=1.8,col="#444444")
    box(xs[i], vy, vw, 8, GREY, GREY_E, labs[[i]][1], labs[[i]][2], cexh = 0.8, cexs = 0.72)
  }
  ## classification
  box(50, 20, 48, 7.5, GREEN, GREEN_E, "Classify against the GTEx baseline", "baseline vs disease-emergent")
  for (i in seq_along(xs)) arrows(xs[i], vy-4, 50, 24, length = 0.06, lwd = 1.3, col = "#666666")
  ## axes
  box(28, 7, 42, 9, BLUE,  BLUE_E,  "Male MAIT axis", "receptor SLC4A10/TRAV1-2; baseline", cexs = 0.76)
  box(74, 7, 42, 9, ORANGE, ORANGE_E, "Female Treg / cDC / CD8 axis", "regulatory/APC; disease-emergent", cexs = 0.76)
  arr(46,16.5,32,11.5); arr(54,16.5,70,11.5)
}

save_one <- function(dev_open, ext) {
  f <- file.path(outdir, paste0("SuppFig_S1_workflow.", ext)); dev_open(f); draw(); dev.off(); message("wrote ", f)
}
save_one(function(f) pdf(f, width = 9.6, height = 9.6), "pdf")
save_one(function(f) {
  if (capabilities("tiff")) tiff(f, width = 1920, height = 1920, res = 200, compression = "lzw") else png(sub("tiff$","png",f), 1920,1920,res=200)
}, "tiff")
