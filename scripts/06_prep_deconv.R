## 06_prep_deconv.R -----------------------------------------------------------
## Prepare per-cohort SYMBOL-level, LINEAR-scale matrices for immunedeconv.
## Maps Entrez/Ensembl -> HGNC symbol (org.Hs.eg.db), CPM-normalises counts,
## de-logs the microarray, collapses duplicate symbols. Writes expr_sym_<GSE>.tsv.
## ---------------------------------------------------------------------------
source("00_config.R")
suppressMessages({library(AnnotationDbi); library(org.Hs.eg.db); library(data.table)})

## per-cohort: gene-id type + scale
CFG <- list(
  GSE130970 = list(idtype = "entrez",  scale = "counts"),
  GSE162694 = list(idtype = "ensembl", scale = "counts"),
  GSE89632  = list(idtype = "symbol",  scale = "log2array"),
  GSE135251 = list(idtype = "ensembl", scale = "counts")
)

to_symbol <- function(ids, idtype) {
  stripped <- sub("\\.\\d+$", "", ids)                 # strip Ensembl version
  if (idtype == "symbol") return(setNames(toupper(ids), ids))
  key <- if (idtype == "entrez") "ENTREZID" else "ENSEMBL"
  m <- suppressMessages(mapIds(org.Hs.eg.db, unique(stripped), "SYMBOL", key, multiVals = "first"))
  out <- toupper(m[stripped]); names(out) <- ids
  out[!is.na(out)]
}

for (gid in names(CFG)) {
  idtype <- CFG[[gid]]$idtype; scale <- CFG[[gid]]$scale
  fp <- file.path(REALDIR, paste0("expr_", gid, ".csv"))
  if (!file.exists(fp)) { cat("  missing", fp, "\n"); next }
  M <- as.data.frame(fread(fp)); rn <- as.character(M[[1]]); M[[1]] <- NULL
  M <- as.matrix(M); rownames(M) <- rn; storage.mode(M) <- "double"
  M <- M[rowSums(!is.na(M)) > 0, , drop = FALSE]
  cat(sprintf("== %s == raw %dx%d idtype=%s scale=%s\n", gid, nrow(M), ncol(M), idtype, scale))

  ## scale -> linear
  if (scale == "counts") {
    cs <- colSums(M, na.rm = TRUE); cs[cs == 0] <- NA
    M <- sweep(M, 2, cs, "/") * 1e6                     # CPM (linear)
    collapse <- "sum"
  } else if (scale == "log2array") {
    M <- 2^M                                            # de-log to linear intensity
    collapse <- "mean"
  } else collapse <- "sum"

  ## map to symbols
  smap <- to_symbol(rownames(M), idtype)
  M <- M[rownames(M) %in% names(smap), , drop = FALSE]
  syms <- smap[rownames(M)]
  ## collapse duplicate symbols
  agg <- if (collapse == "sum") rowsum(M, syms) else {
    cnt <- as.vector(table(syms)[rownames(rowsum(M, syms))]); rowsum(M, syms) / cnt
  }
  agg <- agg[!rownames(agg) %in% c("", "NA", "NONE", "NAN"), , drop = FALSE]
  out <- file.path(REALDIR, paste0("expr_sym_", gid, ".tsv"))
  fwrite(data.frame(symbol = rownames(agg), agg, check.names = FALSE), out, sep = "\t")
  cat(sprintf("   -> %s  symbols=%d  samples=%d  e.g. %s  colmax=%.1f\n",
      basename(out), nrow(agg), ncol(agg), paste(head(rownames(agg),4), collapse=","),
      max(agg, na.rm = TRUE)))
}
cat("DONE – symbol-level linear matrices ready for 07_run_deconv.R\n")
