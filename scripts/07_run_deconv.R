## 07_run_deconv.R -----------------------------------------------------------
## Benchmarked deconvolution (xCell + MCP-counter) of the 4 MASLD cohorts, on
## the symbol-level linear matrices from 06_prep_deconv.R. Writes deconv_xcell.csv
## and deconv_mcp.csv, which 08_concordance.R then tests for sex-direction
## agreement with the singscore results.
## ---------------------------------------------------------------------------
source("00_config.R")
suppressMessages({library(xCell); library(data.table)})
have_mcp <- requireNamespace("MCPcounter", quietly = TRUE)
if (have_mcp) suppressMessages(library(MCPcounter))

xl <- list(); ml <- list()
for (g in COHORT_IDS) {
  fp <- file.path(REALDIR, paste0("expr_sym_", g, ".tsv"))
  if (!file.exists(fp)) { cat("missing", fp, "\n"); next }
  M <- as.data.frame(fread(fp)); rn <- as.character(M[[1]]); M[[1]] <- NULL
  M <- as.matrix(M); rownames(M) <- rn; storage.mode(M) <- "double"
  M <- M[!is.na(rownames(M)) & rownames(M) != "", , drop = FALSE]
  isarray <- (g == "GSE89632")

  xc <- tryCatch(xCellAnalysis(M, rnaseq = !isarray),
                 error = function(e) { cat("xcell err", g, conditionMessage(e), "\n"); NULL })
  if (!is.null(xc)) { d <- as.data.frame(t(xc)); d$sample <- rownames(d); d$cohort <- g; xl[[g]] <- d }

  if (have_mcp) {
    mc <- tryCatch(MCPcounter.estimate(M, featuresType = "HUGO_symbols"),
                   error = function(e) { cat("mcp err", g, conditionMessage(e), "\n"); NULL })
    if (!is.null(mc)) { d <- as.data.frame(t(mc)); d$sample <- rownames(d); d$cohort <- g; ml[[g]] <- d }
  }
  cat(g, "done (samples:", ncol(M), ")\n"); flush.console()
}
if (length(xl)) { fwrite(rbindlist(xl, fill = TRUE), file.path(REALDIR, "deconv_xcell.csv"))
  cat("wrote deconv_xcell.csv\n") }
if (length(ml)) { fwrite(rbindlist(ml, fill = TRUE), file.path(REALDIR, "deconv_mcp.csv"))
  cat("wrote deconv_mcp.csv\n") }
cat("DONE\n")
