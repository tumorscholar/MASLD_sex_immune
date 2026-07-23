## sc_build_owncohort.R -------------------------------------------------------
## Provenance for the in-house CITE-seq validation (Fig 7, section 3.7).
## Regenerates data/own_cohort_percell_fractions.csv, the per-patient table that
## Fig 7 is drawn from, so the whole chain (raw object -> fractions -> figure) is
## reproducible from one place.
##
## Two routes are provided:
##   ROUTE A (canonical)  read your full annotated CITE-seq Seurat object and use
##                        the cell-type labels you already curated (RNA + ADT).
##                        This is how the shipped CSV was made and is the route to
##                        use for the paper.
##   ROUTE B (cross-check) gate the four populations from the exported marker
##                        panel with the SAME rule used for the public cohorts
##                        (sc_validate_*.R). No big object needed. Numbers differ
##                        a little from route A because route A also uses protein
##                        (ADT); route B is only a self-contained sanity check.
##
## Output columns (proportions, 0-1):
##   pid, sex, n_cells, n_immune, n_T,
##   MAIT_ofT, Treg_ofT, CD8_ofT, cDC_ofImm, NK_ofImm
## Fig 7 multiplies the fraction columns by 100 for display.
## ===========================================================================

ROUTE <- Sys.getenv("OWNCOHORT_ROUTE", "A")   # "A" = Seurat object, "B" = panel

## ---- shared config ---------------------------------------------------------
REALDIR <- Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta")
OUT_CSV <- file.path(REALDIR, "single_cell", "own_citeseq",
                     "own_cohort_percell_fractions.csv")
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)

## the four headline populations, gated identically to the public-cohort
## validations (sc_validate_andrews.R etc.) so the in-house and public single-cell analyses share the same gating.
MARK <- list(
  Tcell  = c("CD3D","CD3E","CD3G"),
  MAIT   = c("SLC4A10","TRAV1-2"),
  CD8    = c("CD8A","CD8B"),
  Treg   = c("FOXP3"),
  cDC    = c("CD1C","CLEC9A","CLEC10A","FCER1A"),
  pDC    = c("LILRA4","IL3RA","GZMB"),
  NK     = c("NCAM1","KLRF1","NKG7","KLRD1"),
  immune = c("PTPRC"))

## per-patient fractions from per-cell logical gates + a patient/sex vector
summarise_fractions <- function(pid, sex, is_T, is_MAIT, is_CD8, is_Treg,
                                is_cDC, is_NK, is_imm) {
  pids  <- sort(unique(pid))
  frac  <- function(numer, denom) {
    n <- tapply(numer, pid, sum)[pids]; d <- tapply(denom, pid, sum)[pids]
    ifelse(d > 0, n / d, NA)
  }
  data.frame(
    pid       = pids,
    sex       = tapply(sex, pid, function(x) x[1])[pids],
    n_cells   = as.integer(table(pid)[pids]),
    n_immune  = as.integer(tapply(is_imm, pid, sum)[pids]),
    n_T       = as.integer(tapply(is_T,   pid, sum)[pids]),
    MAIT_ofT  = frac(is_MAIT, is_T),
    Treg_ofT  = frac(is_Treg, is_T),
    CD8_ofT   = frac(is_CD8,  is_T),
    cDC_ofImm = frac(is_cDC,  is_imm),
    NK_ofImm  = frac(is_NK,   is_imm),
    row.names = NULL, stringsAsFactors = FALSE)
}

## ===========================================================================
## ROUTE A - full annotated CITE-seq Seurat object (the canonical route)
## ===========================================================================
if (ROUTE == "A") {
  ## ---- EDIT THESE to match your object -----------------------------------
  OBJ_PATH   <- Sys.getenv("OWN_SEURAT",
                  file.path(REALDIR, "single_cell", "own_citeseq", "liver_citeseq.rds"))
  PID_COL    <- "Sample"        # per-patient id column in meta.data
  SEX_COL    <- "sex"           # recorded/curated sex column ("M"/"F")
  CELLTYPE   <- "cell_type"     # your curated cell-type label column
  ## label strings AS THEY APPEAR in CELLTYPE (regex, case-insensitive):
  LAB <- list(
    Tcell = "T cell|CD4|CD8|Treg|MAIT|Tconv|gdT",
    MAIT  = "MAIT",
    CD8   = "CD8",
    Treg  = "Treg|regulatory",
    cDC   = "cDC|conventional dendritic|DC1|DC2",
    NK    = "\\bNK\\b|natural killer",
    immune= "T cell|CD4|CD8|Treg|MAIT|NK|B cell|Plasma|Mono|Macro|Kupffer|DC|dendritic|Mast|Neutro")
  ## -------------------------------------------------------------------------
  suppressMessages({library(Seurat)})
  message("ROUTE A: loading ", OBJ_PATH)
  obj <- readRDS(OBJ_PATH)
  md  <- obj@meta.data
  for (cc in c(PID_COL, SEX_COL, CELLTYPE))
    if (!cc %in% colnames(md)) stop("column not found in meta.data: ", cc)

  pid <- as.character(md[[PID_COL]])
  pid <- sub("^GC-WL-", "", sub("-LIVER$", "", pid))           # 10113-1 etc.
  sex <- toupper(substr(as.character(md[[SEX_COL]]), 1, 1))     # M / F
  ct  <- as.character(md[[CELLTYPE]])
  hit <- function(key) grepl(LAB[[key]], ct, ignore.case = TRUE)

  is_T    <- hit("Tcell")
  is_MAIT <- is_T & hit("MAIT")
  is_CD8  <- is_T & hit("CD8")
  is_Treg <- is_T & hit("Treg")
  is_cDC  <- hit("cDC")
  is_NK   <- hit("NK")
  is_imm  <- hit("immune")

  df <- summarise_fractions(pid, sex, is_T, is_MAIT, is_CD8, is_Treg,
                            is_cDC, is_NK, is_imm)

## ===========================================================================
## ROUTE B - marker-panel cross-check (self-contained, no big object)
## ===========================================================================
} else {
  PANEL   <- Sys.getenv("OWN_PANEL_CSV",
               file.path(REALDIR, "single_cell", "own_citeseq", "liver_percell_rna.csv"))
  SEXMAP  <- Sys.getenv("OWN_SEXMAP_CSV",
               file.path(dirname(OUT_CSV), "own_cohort_sex.csv"))  # pid,sex
  message("ROUTE B: gating marker panel ", PANEL)
  d <- read.csv(PANEL, stringsAsFactors = FALSE, check.names = FALSE)
  names(d)[1] <- "barcode"
  d$pid <- sub("^GC-WL-", "", sub("-LIVER$", "", d$Sample))
  sm  <- read.csv(SEXMAP, stringsAsFactors = FALSE)
  d$sex <- sm$sex[match(d$pid, sm$pid)]
  d <- d[!is.na(d$sex), , drop = FALSE]

  anyMark <- function(gs) {
    gs <- gs[gs %in% names(d)]
    if (!length(gs)) return(rep(FALSE, nrow(d)))
    Reduce(`|`, lapply(gs, function(g) d[[g]] > 0))
  }
  is_T    <- anyMark(MARK$Tcell)
  is_MAIT <- is_T & anyMark(MARK$MAIT)
  is_CD8  <- is_T & anyMark(MARK$CD8)
  is_Treg <- is_T & anyMark(MARK$Treg)
  is_cDC  <- anyMark(MARK$cDC) & !anyMark(MARK$pDC)
  is_NK   <- !is_T & anyMark(MARK$NK)
  is_imm  <- anyMark(MARK$immune) | is_T

  df <- summarise_fractions(d$pid, d$sex, is_T, is_MAIT, is_CD8, is_Treg,
                            is_cDC, is_NK, is_imm)
}

## ---- write + quick read-out ------------------------------------------------
write.csv(df, OUT_CSV, row.names = FALSE)
cat("\nWrote", OUT_CSV, "\n")
cat(sprintf("Patients: %d  (M=%d, F=%d)\n",
            nrow(df), sum(df$sex == "M"), sum(df$sex == "F")))

test1 <- function(v, dir) {
  m <- df[[v]][df$sex == "M"]; f <- df[[v]][df$sex == "F"]
  p <- tryCatch(wilcox.test(m, f)$p.value, error = function(e) NA)
  hi <- if (median(m, na.rm = TRUE) > median(f, na.rm = TRUE)) "M" else "F"
  cat(sprintf("  %-10s median M=%.3f F=%.3f  higher=%s (predicted %s)  p=%.3f\n",
              v, median(m, na.rm = TRUE), median(f, na.rm = TRUE), hi, dir, p))
}
cat("Headline directions (predicted: MAIT male, Treg/CD8/cDC female):\n")
test1("MAIT_ofT",  "M"); test1("Treg_ofT", "F")
test1("CD8_ofT",   "F"); test1("cDC_ofImm","F")
cat("\nNext: point Fig 7 at this CSV via OWN_CITESEQ_CSV, then run 10_figures.R\n")
