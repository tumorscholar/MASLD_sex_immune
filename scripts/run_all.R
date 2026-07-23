## run_all.R ------------------------------------------------------------------
## End-to-end runner for the sex x MASLD hepatic-immune meta-analysis.
## Runs the whole thing in dependency order: data analysis -> tables -> figures.
## Each step prints a banner and is wrapped so one failure does not kill the run;
## the summary at the end shows what passed, failed or was skipped.
##
## Usage (from the scripts/ folder, in RStudio on the HPC or plain Rscript):
##   setwd("scripts"); source("run_all.R")
## or:
##   Rscript run_all.R
##
## Control what runs with environment variables (all default sensibly):
##   MASLD_REALDIR   output/analysis folder   (see 00_config.R)
##   MASLD_SCRATCH   scratch for downloads     (see 00_config.R)
##   RUN_INSTALL=1   run install_packages.R first (off by default)
##   RUN_BULK=1      bulk pipeline 01-09        (on by default)
##   RUN_SC=0        single-cell validations    (off by default: needs big raw data)
##   RUN_TABLES=1    make_tables.R              (on by default)
##   RUN_FIGURES=1   figures 10/11 + dot plots  (on by default)
## ===========================================================================

## make sure we are in the scripts folder (works from Rscript or source())
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA)
if (!is.na(this_file)) setwd(dirname(this_file))
if (!file.exists("00_config.R"))
  stop("run_all.R must be run from the scripts/ folder (00_config.R not found here).")

flag <- function(v, default) { x <- Sys.getenv(v, NA); if (is.na(x)) default else x %in% c("1","TRUE","true","yes") }
RUN_INSTALL <- flag("RUN_INSTALL", FALSE)
RUN_BULK    <- flag("RUN_BULK",    TRUE)
RUN_SC      <- flag("RUN_SC",      FALSE)
RUN_TABLES  <- flag("RUN_TABLES",  TRUE)
RUN_FIGURES <- flag("RUN_FIGURES", TRUE)

results <- list()
step <- function(label, file, cond = TRUE) {
  if (!cond)             { cat(sprintf("\n---- SKIP  %s (disabled)\n", label)); results[[label]] <<- "skip"; return(invisible()) }
  if (!file.exists(file)){ cat(sprintf("\n---- SKIP  %s (%s not found)\n", label, file)); results[[label]] <<- "skip"; return(invisible()) }
  cat(sprintf("\n==== RUN   %s   [%s]\n%s\n", label, file, strrep("-", 60)))
  t0 <- Sys.time()
  ok <- tryCatch({ source(file, local = new.env()); TRUE },
                 error = function(e) { cat("   ERROR:", conditionMessage(e), "\n"); FALSE })
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("---- %s  %s   (%ss)\n", if (ok) "DONE " else "FAIL ", label, dt))
  results[[label]] <<- if (ok) "ok" else "fail"
}

cat("############################################################\n")
cat("#  sex x MASLD hepatic-immune meta-analysis - full run\n")
cat("#  bulk=", RUN_BULK, " single-cell=", RUN_SC, " tables=", RUN_TABLES,
    " figures=", RUN_FIGURES, "\n", sep = "")
cat("############################################################\n")

## ---- 0. setup --------------------------------------------------------------
step("install packages", "install_packages.R", RUN_INSTALL)

## ---- 1. bulk pipeline (public MASLD cohorts) -------------------------------
step("01 sex assignment",      "01_sex_assignment.R", RUN_BULK)
step("02 build matrix",        "02_build_matrix.R",   RUN_BULK)
step("03 sex main effect",     "03_maineffect.R",     RUN_BULK)
step("04 GTEx control",        "04_gtex_control.R",   RUN_BULK)
step("05 export expression",   "05_export_expr.R",    RUN_BULK)
step("06 prep deconvolution",  "06_prep_deconv.R",    RUN_BULK)
step("07 run deconvolution",   "07_run_deconv.R",     RUN_BULK)
step("08 concordance",         "08_concordance.R",    RUN_BULK)
step("09 deconfound",          "09_deconfound.R",     RUN_BULK)

## ---- 2. single-cell validations (need downloaded raw data) -----------------
##  each build_* makes a Seurat object; each validate_* gates + tests it.
step("sc build Andrews",       "sc_build_andrews.R",       RUN_SC)
step("sc validate Andrews",    "sc_validate_andrews.R",    RUN_SC)
step("sc build Guilliams",     "sc_build_guilliams.R",     RUN_SC)
step("sc validate Guilliams",  "sc_validate_guilliams.R",  RUN_SC)
step("sc Guilliams clean",     "sc_guilliams_clean.R",     RUN_SC)
step("sc build Ramachandran",  "sc_build_ramachandran.R",  RUN_SC)
step("sc validate Ramachandran","sc_validate_ramachandran.R", RUN_SC)
step("sc validate HLiCA",      "sc_validate_hlica.R",      RUN_SC)
step("sc meta forest",         "sc_meta_forest.R",         RUN_SC)
step("sc MAIT functional",     "sc_mait_functional.R",     RUN_SC)
step("sc build in-house (Fig 7 provenance)", "sc_build_owncohort.R", RUN_SC)

## ---- 3. tables -------------------------------------------------------------
step("make tables",            "make_tables.R",       RUN_TABLES)

## ---- 4. figures ------------------------------------------------------------
step("10 main figures (Fig 1-7)", "10_figures.R",     RUN_FIGURES)
step("Fig 7 combined single-cell dot plot","sc_meta_dotplot.R", RUN_FIGURES)
step("11 supplementary figures",  "11_supp_figures.R", RUN_FIGURES)

## ---- summary ---------------------------------------------------------------
cat("\n############################################################\n")
cat("#  RUN SUMMARY\n")
cat("############################################################\n")
for (nm in names(results))
  cat(sprintf("  %-4s  %s\n", toupper(results[[nm]]), nm))
nfail <- sum(unlist(results) == "fail")
cat(sprintf("\n%d ok, %d failed, %d skipped\n",
    sum(unlist(results)=="ok"), nfail, sum(unlist(results)=="skip")))
if (nfail > 0) cat("Some steps failed - check the ERROR lines above.",
                   "Usually a missing input (raw data not downloaded) or a package to install.\n")
