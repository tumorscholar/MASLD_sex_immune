## cohort_scout.R -------------------------------------------------------------
## Scouts candidate public MASLD/NAFLD liver bulk cohorts for inclusion in the
## meta-analysis. For each GEO series it reports: sample count, platform, whether
## the sex-chromosome genes are measurable (so expression-based sex assignment can
## work), and it DUMPS the sample-characteristic fields so we can see how sex,
## disease and fibrosis stage are encoded. Run this, paste the output back, and the
## cohorts that (a) load, (b) carry sex genes and (c) have a codeable fibrosis /
## disease field get wired into 00_config.R + geo_loaders.R.
##
## This is a read-only diagnostic; it does not change the pipeline.
## ===========================================================================
source("00_config.R")
suppressMessages(library(GEOquery))

## edit this list freely; these are the initial candidates
CANDIDATES <- c("GSE48452",   # Ahrens 2013, array, control/steatosis/NASH
                "GSE49541",   # Moylan 2014, array, mild vs advanced fibrosis
                "GSE126848",  # Suppli 2019, RNA-seq, normal/obese/NAFLD/NASH
                "GSE167523",  # NAFLD RNA-seq, steatosis vs NASH
                "GSE83452")   # Lefebvre 2017, array, NASH (note: paired samples)
SEXG <- c("XIST","RPS4Y1","DDX3Y","KDM5D","UTY","EIF1AY")

for (g in CANDIDATES) {
  cat("\n##################### ", g, " #####################\n")
  e <- tryCatch(getGEO(g, destdir = GEO_CACHE, GSEMatrix = TRUE, getGPL = FALSE),
                error = function(err) { cat("  getGEO failed:", conditionMessage(err), "\n"); NULL })
  if (is.null(e)) next
  eset <- e[[1]]
  ph <- Biobase::pData(eset)
  cat("  samples:", nrow(ph), "| platform:", Biobase::annotation(eset), "\n")

  ## expression + sex-gene check (arrays carry expression here; RNA-seq usually don't)
  X <- tryCatch(Biobase::exprs(eset), error = function(err) matrix(nrow = 0, ncol = 0))
  if (nrow(X) > 0) {
    fx <- Biobase::fData(eset)
    symcol <- names(fx)[tolower(names(fx)) %in% c("gene symbol","gene_symbol","symbol","genesymbol")][1]
    sym <- if (!is.na(symcol)) toupper(as.character(fx[[symcol]])) else toupper(rownames(X))
    present <- SEXG[SEXG %in% sym]
    cat("  expression:", nrow(X), "x", ncol(X), "| sex genes present:",
        if (length(present)) paste(present, collapse = ",") else "NONE (check platform annotation)", "\n")
  } else {
    cat("  expression: not in the series matrix (likely RNA-seq; counts are in the",
        "supplementary files and would load via geo_loaders' rnaseq path)\n")
  }

  ## dump the parsed sample-characteristic fields (how sex / disease / fibrosis are encoded)
  cc <- names(ph)[grepl(":ch1$", names(ph)) | grepl("^characteristics", names(ph), ignore.case = TRUE)]
  if (!length(cc)) cc <- names(ph)[seq_len(min(6, ncol(ph)))]
  for (col in unique(cc)) {
    u <- unique(as.character(ph[[col]]))
    cat(sprintf("  [%s]  (%d unique)  %s\n", col, length(u), paste(head(u, 8), collapse = " | ")))
  }
}
cat("\nDone. Paste the whole output back. Cohorts that load, carry sex genes, and have",
    "\na codeable fibrosis/disease field will be added to the pipeline.\n")
