## make_tables.R --------------------------------------------------------------
## Builds the manuscript tables from the pipeline output CSVs, so every number
## in the paper traces back to a script. Writes each table as a CSV under
## $MASLD_REALDIR/tables/ and, if openxlsx is installed, a single workbook
## MASLD_sex_tables.xlsx with one sheet per table.
##
##   Table 1   cohort / dataset characteristics   (analysis_matrix.csv)
##   Table 2   sex main effect per readout        (maineffect_results.csv)
##   Table 3   disease-specificity vs GTEx        (gtex_healthy_sex.csv)
##   Table 4   deconvolution concordance          (deconv_concordance.csv)
##   Table 5   single-cell validation summary     (per-cohort + in-house fractions)
##   Table 5b  single-cell cohort donor counts    (N, M/F per cohort; for §3.7 / Fig 7)
##
## Run after 03 (main effect), 04 (GTEx), 08 (concordance) and the single-cell
## validations.
## ===========================================================================
RD  <- Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta")
SC  <- file.path(RD, "single_cell")
TAB <- file.path(RD, "tables"); dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
rd  <- function(p) if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else { cat("missing:", p, "\n"); NULL }

LAB <- c(ct_MAIT="MAIT", ct_MAITspec="MAIT (SLC4A10/TRAV1-2)", ct_MAITpromisc="MAIT (promiscuous)",
  ct_Treg="Treg", ct_CD8T="CD8 T", ct_DC="Dendritic (cDC)", ct_Tcell="T cell", ct_NK="NK",
  ct_Bcell="B cell", ct_Plasma="Plasma", ct_MonoMac="Mono/Mac", ct_Neutrophil="Neutrophil",
  st_Th1="Th1", st_Th17="Th17", st_CD8_cytotox="CD8 cytotoxic", st_cytotoxCD4="Cytotoxic CD4",
  st_senescence="Senescence", st_exhaustion="Exhaustion", st_Tpex="Tpex", st_Trm="Trm")
lab  <- function(r) ifelse(r %in% names(LAB), LAB[r], r)
star <- function(f) vapply(f, function(x) if (is.na(x)) "" else if (x<1e-3) "***"
                           else if (x<1e-2) "**" else if (x<0.05) "*" else "ns", character(1))
tables <- list()

## ---- Table 1  cohort characteristics ---------------------------------------
mtx <- rd(file.path(RD, "analysis_matrix.csv"))
if (!is.null(mtx)) {
  mtx$fib <- if ("fibrosis_stage" %in% names(mtx)) mtx$fibrosis_stage else NA
  fo <- function(v) { m <- regmatches(v, regexpr("F([0-4])", v))
                      ifelse(nchar(m) > 0, as.integer(sub("F","",m)), NA) }
  mtx$fibn <- suppressWarnings(fo(as.character(mtx$fib)))
  t1 <- do.call(rbind, lapply(split(mtx, mtx$cohort), function(d) data.frame(
    Cohort      = d$cohort[1],
    Platform    = if ("platform" %in% names(d)) d$platform[1] else NA,
    N           = nrow(d),
    Male        = sum(d$sex_assigned == "M", na.rm = TRUE),
    Female      = sum(d$sex_assigned == "F", na.rm = TRUE),
    Pct_male    = round(100 * mean(d$sex_assigned == "M", na.rm = TRUE), 1),
    F0_F1       = sum(d$fibn %in% 0:1, na.rm = TRUE),
    F2          = sum(d$fibn == 2, na.rm = TRUE),
    F3_F4       = sum(d$fibn %in% 3:4, na.rm = TRUE),
    row.names = NULL)))
  t1 <- rbind(t1, data.frame(Cohort="TOTAL", Platform="",
    N=sum(t1$N), Male=sum(t1$Male), Female=sum(t1$Female),
    Pct_male=round(100*sum(t1$Male)/sum(t1$Male+t1$Female),1),
    F0_F1=sum(t1$F0_F1), F2=sum(t1$F2), F3_F4=sum(t1$F3_F4)))
  tables[["Table1_cohorts"]] <- t1
}

## ---- Table 2  sex main effect (fibrosis-adjusted, ordinal coding) -----------
me <- rd(file.path(RD, "maineffect_results.csv"))
cc <- rd(file.path(RD, "maineffect_concordance.csv"))
if (!is.null(me)) {
  m <- me[me$coding == "fib_ord", ]
  head_set <- if (!is.null(cc)) cc$readout[cc$headline %in% c(TRUE,"True","TRUE")] else character(0)
  t2 <- data.frame(
    Readout   = lab(m$readout),
    N         = m$n,
    Beta_sexM = round(m$beta_sexM, 3),
    SE        = round(m$se, 3),
    Direction = ifelse(m$beta_sexM > 0, "male", "female"),
    P         = signif(m$p, 3),
    FDR       = signif(m$fdr, 3),
    Sig       = star(m$fdr),
    LOCO_consistent = m$loco_sign_consistent,
    LOCO_n_sig = m$loco_n_sig,
    Headline  = m$readout %in% head_set,
    row.names = NULL)
  t2 <- t2[order(t2$FDR), ]
  tables[["Table2_main_effect"]] <- t2
}

## ---- Table 3  disease-specificity (MASLD vs disease-free GTEx) --------------
gt <- rd(file.path(RD, "gtex_healthy_sex.csv"))
if (is.null(gt)) gt <- rd(file.path(dirname(RD), "results_R/gtex_healthy_sex.csv"))
if (is.null(gt)) gt <- rd("../results_R/gtex_healthy_sex.csv")
if (!is.null(gt)) {
  t3 <- data.frame(
    Readout            = lab(gt$readout),
    MASLD_direction    = gt$MASLD_dir,
    GTEx_beta_adj      = round(gt$beta_cov, 3),
    GTEx_p_adj         = signif(gt$p_cov, 3),
    GTEx_direction     = gt$dir_cov,
    Balanced_consistency = gt$bal_consistency,
    Verdict            = gt$verdict,
    row.names = NULL)
  tables[["Table3_disease_specificity"]] <- t3
}

## ---- Table 4  benchmarked deconvolution concordance ------------------------
dc <- rd(file.path(RD, "deconv_concordance.csv"))
if (!is.null(dc)) {
  keep <- intersect(c("method","deconv_celltype","singscore_beta","deconv_beta"), names(dc))
  if (length(keep) >= 3) {
    t4 <- data.frame(
      Method        = dc$method,
      Cell_type     = dc$deconv_celltype,
      singscore_beta = round(dc$singscore_beta, 3),
      deconv_beta    = round(dc$deconv_beta, 3),
      Concordant     = ifelse(sign(dc$singscore_beta) == sign(dc$deconv_beta), "yes", "no"),
      row.names = NULL)
    tables[["Table4_deconvolution_concordance"]] <- t4
  }
}

## ---- Table 5  single-cell validation summary -------------------------------
## one row per cohort x readout: median M, median F, direction, predicted, p.
EXPECT <- c(MAIT="M", Treg="F", CD8="F", cDC="F")
sc_rows <- list()

add_cohort <- function(cohort_label, df, cols) {
  df <- df[df$sex %in% c("M","F"), , drop = FALSE]
  if (sum(df$sex=="M") < 2 || sum(df$sex=="F") < 2) return(invisible())
  for (nm in names(cols)) {
    col <- cols[[nm]]; if (!col %in% names(df)) next
    v <- df[[col]]; if (max(v, na.rm=TRUE) <= 1.5) v <- v*100   # proportion -> %
    m <- v[df$sex=="M"]; f <- v[df$sex=="F"]
    p <- tryCatch(wilcox.test(m,f)$p.value, error=function(e) NA)
    hi <- if (median(m,na.rm=TRUE) > median(f,na.rm=TRUE)) "M" else "F"
    sc_rows[[length(sc_rows)+1]] <<- data.frame(
      Cohort=cohort_label, Readout=nm, N_M=sum(df$sex=="M"), N_F=sum(df$sex=="F"),
      Median_M=round(median(m,na.rm=TRUE),2), Median_F=round(median(f,na.rm=TRUE),2),
      Higher_in=hi, Predicted=EXPECT[[nm]], Matches=ifelse(hi==EXPECT[[nm]],"yes","no"),
      P=signif(p,3), row.names=NULL)
  }
}

## in-house cohort (ships in data/)
own_paths <- c(file.path(SC,"own_citeseq/own_cohort_percell_fractions.csv"),
               file.path(dirname(RD),"data/own_cohort_percell_fractions.csv"),
               "../data/own_cohort_percell_fractions.csv",
               "data/own_cohort_percell_fractions.csv")
own_p <- own_paths[file.exists(own_paths)][1]
if (length(own_p) && !is.na(own_p)) add_cohort("In-house CITE-seq", read.csv(own_p, stringsAsFactors=FALSE),
  list(MAIT="MAIT_ofT", Treg="Treg_ofT", CD8="CD8_ofT", cDC="cDC_ofImm"))

## public cohorts (per-cohort validation fraction CSVs, if present)
pub <- list(
  list(label="Andrews (healthy)",     csv="andrews2024/out/Andrews2024_per_donor_fractions.csv",   sub=c(group="healthy")),
  list(label="HLiCA (healthy)",       csv="hlica/out/HLiCA_per_patient.csv",                        sub=NULL),
  list(label="Guilliams (lean)",      csv="guilliams/out/Guilliams_clean_per_patient.csv",          sub=c(diet="Lean")),
  list(label="Ramachandran (cirrh.)", csv="ramachandran/out/Ramachandran2019_per_donor_fractions.csv", sub=c(group="cirrhotic")))
for (co in pub) {
  fp <- file.path(SC, co$csv); if (!file.exists(fp)) { cat("missing:", co$csv, "\n"); next }
  x <- read.csv(fp, stringsAsFactors=FALSE); names(x)[names(x)=="n_Tcells"] <- "n_T"
  if (!is.null(co$sub) && names(co$sub) %in% names(x)) x <- x[as.character(x[[names(co$sub)]])==unname(co$sub),,drop=FALSE]
  if ("n_T" %in% names(x)) x <- x[!is.na(x$n_T) & x$n_T >= 50, ]
  add_cohort(co$label, x, list(MAIT="MAIT_pctT", Treg="Treg_pctT", CD8="CD8_pctT", cDC="cDC_pctImm"))
}
if (length(sc_rows)) {
  tables[["Table5_singlecell_validation"]] <- do.call(rbind, sc_rows)
  ## Table 5b: one row per single-cell cohort with donor counts (for §3.7 / Fig 7)
  scv <- do.call(rbind, sc_rows)
  cc <- unique(scv[, c("Cohort","N_M","N_F")])
  cc$N_total <- cc$N_M + cc$N_F
  tables[["Table5b_singlecell_cohorts"]] <- cc[, c("Cohort","N_total","N_M","N_F")]
}

## ---- write CSVs + optional workbook ----------------------------------------
for (nm in names(tables)) {
  write.csv(tables[[nm]], file.path(TAB, paste0(nm, ".csv")), row.names = FALSE)
  cat("wrote", nm, "(", nrow(tables[[nm]]), "rows )\n")
}
if (requireNamespace("openxlsx", quietly = TRUE) && length(tables)) {
  wb <- openxlsx::createWorkbook()
  for (nm in names(tables)) {
    sh <- substr(nm, 1, 31); openxlsx::addWorksheet(wb, sh)
    openxlsx::writeData(wb, sh, tables[[nm]])
    openxlsx::setColWidths(wb, sh, cols = seq_len(ncol(tables[[nm]])), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, file.path(TAB, "MASLD_sex_tables.xlsx"), overwrite = TRUE)
  cat("wrote MASLD_sex_tables.xlsx (", length(tables), "sheets )\n")
} else cat("openxlsx not installed - wrote CSVs only (install.packages('openxlsx') for the workbook)\n")

cat("\nAll tables written to", TAB, "\n")
