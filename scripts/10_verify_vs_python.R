## 10_verify_vs_python.R ------------------------------------------------------
## Verify the R pipeline reproduces the original Python results.
## Run this AFTER a full, fresh 01->09 R run that wrote into a CLEAN output
## directory (R_DIR). It compares each R-generated CSV against the archived
## Python CSV (PY_DIR) and reports, per readout, the difference in the sex
## effect (beta) and whether the direction (sign) agrees.
##
## Betas should match very closely (the scoring/z-scoring/OLS are identical
## ports). P-values may differ slightly for the mixed-model outputs, because
## R (lmerTest, Satterthwaite d.f.) and Python (statsmodels, Wald) compute them
## differently by design — so BETA agreement and SIGN agreement are the real
## test, and p/FDR are shown for information.
## ---------------------------------------------------------------------------
source("00_config.R")

## the fresh R output directory (set MASLD_REALDIR to this before running 01-09)
R_DIR  <- Sys.getenv("MASLD_REALDIR", REALDIR)
## the original Python outputs to compare against
PY_DIR <- Sys.getenv("MASLD_PYDIR", "/data/Blizard-AlazawiLab/rk/MASLD_sex_meta")

BETA_TOL <- 0.02   # |R beta - Python beta| below this = match

cat("R  outputs:", R_DIR,  "\n")
cat("PY outputs:", PY_DIR, "\n")
if (normalizePath(R_DIR, mustWork = FALSE) == normalizePath(PY_DIR, mustWork = FALSE))
  stop("R_DIR and PY_DIR are the SAME folder — point the R run at a fresh directory first ",
       "(e.g. Sys.setenv(MASLD_REALDIR='/data/home/hdx044/MASLD_sex_immune/results_R')).")

## generic comparator: merge on key cols, report diff on numeric cols + sign agreement
compare_csv <- function(fname, key, num, signcol) {
  fr <- file.path(R_DIR, fname); fp <- file.path(PY_DIR, fname)
  if (!file.exists(fr) || !file.exists(fp)) {
    cat(sprintf("[skip ] %-28s missing in %s\n", fname,
                if (!file.exists(fr)) "R output" else "Python output")); return(NULL)
  }
  a <- read.csv(fr, stringsAsFactors = FALSE, check.names = FALSE)
  b <- read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE)
  key <- key[key %in% names(a) & key %in% names(b)]
  m <- merge(a, b, by = key, suffixes = c(".R", ".PY"))
  if (!nrow(m)) { cat(sprintf("[warn ] %-28s no matching rows on key\n", fname)); return(NULL) }
  worst_beta <- NA; sign_ok <- NA
  for (col in num) {
    cR <- paste0(col, ".R"); cP <- paste0(col, ".PY")
    if (!(cR %in% names(m) && cP %in% names(m))) next
    d <- abs(suppressWarnings(as.numeric(m[[cR]])) - suppressWarnings(as.numeric(m[[cP]])))
    md <- max(d, na.rm = TRUE)
    if (col == signcol) worst_beta <- md
    cat(sprintf("        %-28s %-14s max|Δ| = %.4f\n", fname, col, md))
  }
  if (signcol %in% num) {
    sR <- sign(suppressWarnings(as.numeric(m[[paste0(signcol,".R")]])))
    sP <- sign(suppressWarnings(as.numeric(m[[paste0(signcol,".PY")]])))
    sign_ok <- mean(sR == sP, na.rm = TRUE)
  }
  verdict <- if (!is.na(worst_beta) && worst_beta <= BETA_TOL &&
                 (is.na(sign_ok) || sign_ok == 1)) "PASS" else "CHECK"
  cat(sprintf("  [%s] %-26s n=%d  worst beta Δ=%.4f  sign agreement=%s\n\n",
              verdict, fname, nrow(m), worst_beta %||% NA,
              if (is.na(sign_ok)) "-" else sprintf("%.0f%%", 100*sign_ok)))
  data.frame(file = fname, n = nrow(m), worst_beta_diff = worst_beta,
             sign_agreement = sign_ok, verdict = verdict, stringsAsFactors = FALSE)
}
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

cat("\n==== comparing R vs Python outputs (beta tolerance =", BETA_TOL, ") ====\n\n")
res <- rbind(
  compare_csv("maineffect_results.csv",  c("readout","coding"),
              c("beta_sexM","p","fdr"), "beta_sexM"),
  compare_csv("gtex_healthy_sex.csv",    "readout",
              c("beta","p"), "beta"),
  compare_csv("deconv_concordance.csv",  c("method","singscore_readout","deconv_celltype"),
              c("deconv_beta"), "deconv_beta"),
  compare_csv("deconf_within_cohort.csv","readout",
              c("beta_unadj","beta_adj","p_unadj","p_adj"), "beta_unadj")
)

## also spot-check the singscore matrix itself (the scoring port), keyed on sample
amf_R <- file.path(R_DIR, "analysis_matrix.csv"); amf_P <- file.path(PY_DIR, "analysis_matrix.csv")
if (file.exists(amf_R) && file.exists(amf_P)) {
  A <- read.csv(amf_R, stringsAsFactors = FALSE); B <- read.csv(amf_P, stringsAsFactors = FALSE)
  m <- merge(A, B, by = "gsm", suffixes = c(".R",".PY"))
  scols <- grep("^(ct_|st_)", names(A), value = TRUE)
  worst <- 0
  for (c in scols) {
    cR <- paste0(c,".R"); cP <- paste0(c,".PY")
    if (cR %in% names(m) && cP %in% names(m))
      worst <- max(worst, max(abs(m[[cR]] - m[[cP]]), na.rm = TRUE))
  }
  cat(sprintf("  [%s] analysis_matrix.csv        n=%d samples  worst singscore Δ=%.5f\n",
              if (worst < 1e-3) "PASS" else "CHECK", nrow(m), worst))
}

cat("\n==== SUMMARY ====\n")
if (!is.null(res)) print(res, row.names = FALSE)
cat("\nPASS = R reproduces Python (betas within tolerance, directions identical).\n")
cat("CHECK = investigate; a small p-value gap on the mixed-model files is expected\n")
cat("(lmerTest Satterthwaite vs statsmodels Wald) and is not itself a failure.\n")
