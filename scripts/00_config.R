## 00_config.R --------------------------------------------------------------
## Shared paths, marker/signature definitions, and helper functions for the
## sex x MASLD hepatic-immune meta-analysis. Every other script starts with
##   source("00_config.R")
## Edit ONLY the paths in this file to run on your own machine / HPC.
## ---------------------------------------------------------------------------

## ---- paths ----------------------------------------------------------------
## By default these are placeholder paths; set them to your data and scratch
## On a laptop, point them anywhere writable (they are created if missing).
REALDIR <- Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta")
SCRATCH <- Sys.getenv("MASLD_SCRATCH", "/path/to/scratch/masld_meta")

## FAIL LOUDLY if the paths were never set (the #1 cause of cryptic downstream
## "cannot open connection" errors). Set them before sourcing, ideally once in
## ~/.Renviron:   MASLD_REALDIR=/your/output/dir
##                MASLD_SCRATCH=/your/scratch/dir
if (grepl("^/path/to", REALDIR) || grepl("^/path/to", SCRATCH))
  stop("\n  MASLD_REALDIR / MASLD_SCRATCH are still the PLACEHOLDER paths, so nothing will be found.\n",
       "  Set them before sourcing 00_config.R. Easiest (persists across sessions) - add to ~/.Renviron:\n",
       "      MASLD_REALDIR=/data/Blizard-AlazawiLab/rk/MASLD_sex_meta\n",
       "      MASLD_SCRATCH=/data/Blizard-AlazawiLab/rk/scratch/masld_meta\n",
       "  then restart R. Or per-session:  Sys.setenv(MASLD_REALDIR=..., MASLD_SCRATCH=...)\n",
       "  Current values -> REALDIR=", REALDIR, "  SCRATCH=", SCRATCH, "\n", call. = FALSE)

GEO_CACHE <- file.path(SCRATCH, "geo_cache")
FIGDIR  <- file.path(REALDIR, "figures")
for (d in c(REALDIR, SCRATCH, GEO_CACHE, FIGDIR))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- cohorts --------------------------------------------------------------
COHORTS <- list(
  list(gse = "GSE130970", type = "rnaseq"),
  list(gse = "GSE162694", type = "rnaseq"),
  list(gse = "GSE89632",  type = "array"),
  list(gse = "GSE135251", type = "rnaseq"),
  ## public MASLD/NAFLD liver cohorts (selected after scouting + sex-assignment QC):
  list(gse = "GSE167523", type = "rnaseq"),   # n=98, NAFL vs NASH, sex 0 ambiguous, 94.9% concordant
  list(gse = "GSE48452",  type = "array")     # n=73 (50 confident), fibrosis 0-4 + NAS, 98.0% concordant
  ## wave 2 (need per-cohort handling, not yet wired):
  ##   GSE126848 rnaseq n=57 - count-matrix columns are bare sample codes that do not map to GSM
  ##             accessions (needs a code->GSM map before it will merge); has healthy/obese controls
  ##   GSE49541  array n=72, mild(F0-1)/advanced(F3-4) fibrosis, NO recorded sex (expression-only)
  ##   GSE83452  array n=231 LONGITUDINAL (baseline+follow-up, bariatric/diet) - baseline-only before use
)
COHORT_IDS <- vapply(COHORTS, function(x) x$gse, character(1))

## ---- sex markers ----------------------------------------------------------
FEMALE_MARKERS <- c("XIST")
MALE_MARKERS   <- c("RPS4Y1","DDX3Y","EIF1AY","UTY","KDM5D",
                    "USP9Y","NLGN4Y","ZFY","TXLNGY")
SEX_MARGIN <- 0.5   # |z_Y - z_XIST| below this = "Ambiguous" (dropped)

## ---- cell-type identity signatures (up-only) ------------------------------
CELLTYPE <- list(
  ct_CD8T       = c("CD8A","CD8B"),
  ct_Tcell      = c("CD3D","CD3E","CD3G","TRAC"),
  ct_NK         = c("NCAM1","KLRD1","NKG7","KLRF1","NCR1","GNLY"),
  ct_Bcell      = c("CD19","MS4A1","CD79A","CD79B","BANK1"),
  ct_Plasma     = c("MZB1","XBP1","SDC1","PRDM1","DERL3"),
  ct_MonoMac    = c("CD68","CD14","LYZ","CSF1R","ITGAM","FCGR3A"),
  ct_DC         = c("FLT3","CLEC9A","BATF3","CD1C","CLEC10A"),
  ct_Neutrophil = c("FCGR3B","CSF3R","S100A8","S100A9","CXCR2"),
  ct_Treg       = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"),
  ct_MAIT       = c("SLC4A10","KLRB1","ZBTB16","RORC","TRAV1-2"),
  ## receptor-identity MAIT (invariant TCR: SLC4A10 + TRAV1-2) vs promiscuous
  ## CD161/RORgt MAIT-like cells. The receptor-specific readout is the cleanest
  ## male-biased signal and must be scored here so it reaches analysis_matrix.csv
  ## (the figure scripts, 13_meta_random and make_tables all expect these columns).
  ct_MAITspec    = c("SLC4A10","TRAV1-2"),
  ct_MAITpromisc = c("KLRB1","RORC","ZBTB16")
)

## ---- functional-state signatures (list(up, down)) -------------------------
STATE <- list(
  st_CD8_cytotox = list(up = c("GZMA","GZMB","GZMH","PRF1","GNLY","NKG7","KLRD1",
                               "FGFBP2","KLRG1","CST7","CTSW"), down = character(0)),
  st_exhaustion  = list(up = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1",
                               "TOX","CD160","BTLA","LAYN","TNFRSF9","VSIR"), down = character(0)),
  st_Tpex        = list(up = c("TCF7","SLAMF6","CXCR5","IL7R","ID3","BCL6"), down = character(0)),
  st_Trm         = list(up = c("CD69","ITGAE","CXCR6","ITGA1","RBPJ","ZNF683"),
                        down = c("KLF2","S1PR1","SELL","CCR7")),
  st_Th1         = list(up = c("TBX21","IFNG","CXCR3","IL12RB2","STAT1"), down = character(0)),
  st_Th17        = list(up = c("RORC","IL17A","IL17F","CCR6","IL23R"), down = character(0)),
  st_cytotoxCD4  = list(up = c("GZMB","GZMH","GZMA","PRF1","NKG7","GNLY"), down = character(0)),
  st_senescence  = list(up = c("CDKN1A","CDKN2A","IL6","CXCL8","IGFBP3","IGFBP7",
                               "SERPINE1","IL1B","IL1A","CCL2","MMP3","TNFRSF1B"), down = character(0))
)

ALL_SYMS <- sort(unique(c(
  MALE_MARKERS, FEMALE_MARKERS,
  unlist(CELLTYPE, use.names = FALSE),
  unlist(lapply(STATE, function(s) c(s$up, s$down)), use.names = FALSE)
)))

## ---- helpers --------------------------------------------------------------

## within-vector z-score (population standard deviation, denominator n)
zscore <- function(x) {
  x  <- suppressWarnings(as.numeric(x))
  mu <- mean(x, na.rm = TRUE)
  sdv <- sqrt(mean((x - mu)^2, na.rm = TRUE))   # population SD
  if (is.na(sdv) || sdv == 0) return(x * 0)
  (x - mu) / sdv
}

## Signature scoring uses the Bioconductor singscore package throughout
## (singscore::rankGenes + simpleScore); see 02_build_matrix.R and
## 04_gtex_control.R. The earlier hand-rolled rank scorer has been removed so a
## single, standard method is used for both MASLD and the GTEx control.

## map a fibrosis-stage label to an ordinal 0-4 (NA if not codeable)
code_fibrosis <- function(s) {
  vapply(s, function(v) {
    if (is.na(v)) return(NA_real_)
    m <- regmatches(v, regexpr("F([0-4])", v))
    if (length(m) == 1 && nchar(m) > 0) return(as.numeric(sub("F", "", m)))
    if (grepl("normal|healthy|\\bhc\\b|control", v, ignore.case = TRUE)) return(0)
    ## disease-category severity proxy, used ONLY as a within-cohort nuisance
    ## covariate where a fibrosis stage is not reported (e.g. NAFL vs NASH, or
    ## mild vs advanced). Kept ordinal so each cohort adjusts on its own severity
    ## axis; readouts are z-scored within cohort so scales are not compared across.
    if (grepl("nash|advanced|severe", v, ignore.case = TRUE)) return(2)
    if (grepl("nafl|steatos|\\bmild\\b|obese", v, ignore.case = TRUE)) return(1)
    NA_real_
  }, numeric(1))
}

## log2(x+1) if the matrix looks like linear counts (max > 40), else pass-through
maybe_log2 <- function(M) {
  if (max(M, na.rm = TRUE) > 40) log2(pmax(M, 0) + 1) else M
}

message("00_config.R loaded | REALDIR=", REALDIR)
