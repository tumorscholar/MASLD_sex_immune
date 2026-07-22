## 05_export_expr.R -----------------------------------------------------------
## Export per-cohort gene x GSM expression matrices on the LINEAR scale, for the
## deconvolution branch (xCell / MCP-counter expect non-log expression).
## Does NOT touch analysis_matrix.csv. Writes expr_<GSE>.csv (native gene IDs).
## ---------------------------------------------------------------------------
source("00_config.R")
source("geo_loaders.R")

for (co in COHORTS) {
  gid <- co$gse; typ <- co$type
  cat("======", gid, "======\n")
  res <- try({
    gse_soft <- getGEO(gid, destdir = GEO_CACHE, GSEMatrix = FALSE)
    if (typ == "array") {
      eset <- getGEO(gid, destdir = GEO_CACHE, GSEMatrix = TRUE, getGPL = TRUE)[[1]]
      M <- array_symbol_matrix(eset)
    } else {
      M <- as.matrix(map_cols_to_gsm(gse_soft, rnaseq_matrix(gid)))
    }
    M <- M[!duplicated(rownames(M)), , drop = FALSE]
    M <- M[rowSums(!is.na(M)) > 0, , drop = FALSE]
    ## keep LINEAR scale (do NOT log) — deconvolution methods expect it
    outdf <- data.frame(gene_id = rownames(M), M, check.names = FALSE)
    fp <- file.path(REALDIR, paste0("expr_", gid, ".csv"))
    write.csv(outdf, fp, row.names = FALSE)
    cat(sprintf("  wrote %s  shape=%dx%d  max=%.1f  id_sample=%s\n",
        basename(fp), nrow(M), ncol(M), max(M, na.rm = TRUE),
        paste(head(rownames(M), 5), collapse=",")))
  }, silent = TRUE)
  if (inherits(res, "try-error")) cat("  !! FAILED", gid, ":", conditionMessage(attr(res,"condition")), "\n")
}
cat("DONE\n")
