## diag_newcohorts.R ----------------------------------------------------------
## Read-only follow-up to 01_sex_assignment.R. Answers three questions for the
## candidate cohorts without re-running the whole assignment:
##   (A) GSE167523 / GSE48452 : do the confident sex calls agree with recorded sex?
##   (B) GSE126848            : why did it fail to load (which stage, what error)?
## Paste the whole output back.
## ===========================================================================
source("00_config.R")
source("geo_loaders.R")

## ---- (A) concordance for the already-assigned new cohorts -------------------
csv <- file.path(REALDIR, "sex_assign_out", "per_sample_sex_calls.csv")
if (file.exists(csv)) {
  A <- read.csv(csv, stringsAsFactors = FALSE)
  for (g in c("GSE167523", "GSE48452")) {
    d <- A[A$gse == g, ]
    if (!nrow(d)) { cat("\n==", g, "== not in per_sample file\n"); next }
    kn <- d[!is.na(d$recorded_sex) & d$sex_assigned != "Ambiguous", ]
    cat(sprintf("\n== %s ==  total=%d  confident=%d  ambiguous=%d\n",
        g, nrow(d), nrow(kn), sum(d$sex_assigned == "Ambiguous")))
    if (nrow(kn)) {
      cat(sprintf("  concordance vs recorded sex: %.1f%%\n",
          100 * mean(kn$recorded_sex == kn$sex_assigned)))
      print(table(recorded = kn$recorded_sex, assigned = kn$sex_assigned))
    } else cat("  no recorded sex available to compare\n")
  }
} else cat("per_sample_sex_calls.csv not found at", csv, "\n")

## ---- (B) why GSE126848 failed: probe each loader stage ----------------------
cat("\n########## GSE126848 load diagnosis ##########\n")
ddir <- file.path(GEO_CACHE, "GSE126848")
if (!dir.exists(ddir) || !length(list.files(ddir, recursive = TRUE)))
  try(getGEOSuppFiles("GSE126848", baseDir = GEO_CACHE, makeDirectory = TRUE), silent = TRUE)
ff <- list.files(ddir, recursive = TRUE, full.names = TRUE)
cat("supplementary files:\n"); print(basename(ff))

step <- tryCatch({
  M0 <- rnaseq_matrix("GSE126848")
  cat(sprintf("  [1] rnaseq_matrix OK: %d x %d\n", nrow(M0), ncol(M0)))
  cat("      first row ids : ", paste(head(rownames(M0), 4), collapse = " | "), "\n")
  cat("      first col names: ", paste(head(colnames(M0), 4), collapse = " | "), "\n")
  M1 <- map_cols_to_gsm(getGEO("GSE126848", destdir = GEO_CACHE, GSEMatrix = FALSE), M0)
  cat("  [2] map_cols_to_gsm OK; cols now: ",
      paste(head(colnames(M1), 4), collapse = " | "), "\n")
  ids <- rownames(M1); up <- toupper(ids); ens <- toupper(sub("\\.\\d+$", "", ids))
  ENS <- c(XIST = "ENSG00000229807", RPS4Y1 = "ENSG00000129824",
           DDX3Y = "ENSG00000067048", KDM5D = "ENSG00000012817")
  fsym <- c("XIST","RPS4Y1","DDX3Y","KDM5D")[c("XIST","RPS4Y1","DDX3Y","KDM5D") %in% up]
  fens <- names(ENS)[ENS %in% ens]
  cat("  [3] sex markers present -> by symbol:",
      if (length(fsym)) paste(fsym, collapse = ",") else "none",
      "| by Ensembl:", if (length(fens)) paste(fens, collapse = ",") else "none", "\n")
  "reached end OK (loads fine now)"
}, error = function(e) paste("STOPPED with error:", conditionMessage(e)))
cat("  =>", step, "\n")
cat("\nDone. Paste the whole output back.\n")
