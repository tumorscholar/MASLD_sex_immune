## 04_gtex_control.R ----------------------------------------------------------
## Disease-specificity control. Score the identical immune readouts in GTEx
## disease-free liver (~n=226) and test the sex effect, to ask whether each
## MASLD sex bias is CONSTITUTIONAL (present in healthy liver too) or
## DISEASE-EMERGENT (appears only in MASLD).
##   readout ~ sexM + age   (plain OLS; GTEx is a single cohort)
##
## NOTE (known limitation, flagged by review): this baseline model adjusts only
## for age. GTEx liver is post-mortem, immune-poor and not sex-balanced; the
## planned strengthening is to add technical covariates (RIN = SMRIN,
## ischaemic time = SMTSISCH, centre = SMCENTER) and sex-balance the comparison.
## Those columns live in the SampleAttributes file already downloaded below.
## ---------------------------------------------------------------------------
source("00_config.R")
suppressMessages(library(data.table))

BASE <- "https://storage.googleapis.com/adult-gtex"
FILES <- c(
  attr = paste0(BASE, "/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
  phen = paste0(BASE, "/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"),
  tpm  = paste0(BASE, "/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
)
fetch <- function(k) {
  p <- file.path(SCRATCH, basename(FILES[k]))
  if (!file.exists(p)) {
    message("  downloading ", k)
    download.file(FILES[k], p, mode = "wb", quiet = TRUE)
  }
  p
}

## GTEx-specific readout sets (includes the MAIT receptor-specific / promiscuous splits)
GT_CELLTYPE <- CELLTYPE
GT_CELLTYPE$ct_MAITspec    <- c("SLC4A10","TRAV1-2")
GT_CELLTYPE$ct_MAITpromisc <- c("KLRB1","RORC","ZBTB16")
GT_STATE <- STATE[c("st_CD8_cytotox","st_Th1","st_Th17","st_cytotoxCD4","st_senescence")]
MASLD_DIR <- c(ct_MAIT="MALE", ct_MAITspec="MALE", ct_Treg="female",
               ct_CD8T="female", ct_DC="female", st_Th1="female")

## ---- annotations ----
message("Fetching GTEx annotations...")
attr <- fread(fetch("attr"), sep = "\t")
phen <- fread(fetch("phen"), sep = "\t")
liver <- attr[SMTSD == "Liver", SAMPID]
cat("GTEx liver samples:", length(liver), "\n")
setkey(phen, SUBJID)
sid <- function(s) vapply(strsplit(s, "-"), function(x) paste(x[1:2], collapse = "-"), character(1))

## ---- stream the gene-TPM gct, keeping only liver columns (memory-safe) ----
message("Reading GTEx gene TPM (streaming liver columns only)...")
con <- gzfile(fetch("tpm"), "rt")
readLines(con, 2)                                  # skip the 2 gct header lines
header <- strsplit(readLines(con, 1), "\t")[[1]]
di <- which(header == "Description")
liv <- liver[liver %in% header]
li  <- match(liv, header)
genes <- character(0); mat <- vector("list", 0); k <- 0
repeat {
  chunk <- readLines(con, 2000)
  if (length(chunk) == 0) break
  sp <- strsplit(chunk, "\t")
  genes <- c(genes, vapply(sp, function(p) p[di], character(1)))
  for (p in sp) { k <- k + 1; mat[[k]] <- as.numeric(p[li]) }
}
close(con)
M <- do.call(rbind, mat); rownames(M) <- genes; colnames(M) <- liv
M <- M[!duplicated(rownames(M)), , drop = FALSE]
M <- log2(pmax(M, 0) + 1)
cat("Expression matrix:", nrow(M), "x", ncol(M), "\n")

## ---- score ----
ng <- nrow(M); ranks <- rank_matrix(M)
score <- list()
for (nm in names(GT_CELLTYPE)) {
  s <- singscore_custom(ranks, ng, GT_CELLTYPE[[nm]]); if (!is.null(s)) score[[nm]] <- s
}
for (nm in names(GT_STATE)) {
  s <- singscore_custom(ranks, ng, GT_STATE[[nm]]$up, GT_STATE[[nm]]$down); if (!is.null(s)) score[[nm]] <- s
}
S <- as.data.frame(score)                          # samples x readouts
## singscore_custom returns vectors named by sample (colMeans keeps colnames of M)
S$SAMPID <- rownames(S)

## ---- donor metadata ----
S$subj <- sid(S$SAMPID)
S$SEX  <- phen[S$subj, SEX]
S$AGE  <- phen[S$subj, AGE]
S$sexM <- as.numeric(S$SEX == 1)
agemid <- function(a) vapply(a, function(v) {
  z <- strsplit(v, "-")[[1]]; if (length(z)==2) (as.numeric(z[1])+as.numeric(z[2]))/2 else NA_real_
}, numeric(1))
S$age <- agemid(S$AGE)
S <- S[!is.na(S$sexM), , drop = FALSE]
cat(sprintf("GTEx liver donors: M=%d F=%d\n", sum(S$sexM==1), sum(S$sexM==0)))

## ---- sex effect in disease-free liver ----
cat("\n==== SEX effect in DISEASE-FREE GTEx liver (readout ~ sexM + age) ====\n")
cat("(beta>0 = higher in MALES; compare 'GTEx_dir' vs 'MASLD_dir')\n")
readout_names <- c(names(GT_CELLTYPE), names(GT_STATE))
rows <- list()
for (nm in readout_names) {
  if (!nm %in% names(S)) next
  d <- S[!is.na(S[[nm]]), , drop = FALSE]; d[[nm]] <- zscore(d[[nm]])
  fit <- tryCatch(lm(as.formula(sprintf("%s ~ sexM + age", nm)), d), error = function(e) NULL)
  if (is.null(fit)) next
  co <- summary(fit)$coefficients
  b <- co["sexM","Estimate"]; p <- co["sexM","Pr(>|t|)"]
  gd <- if (b > 0) "MALE" else "female"
  md <- if (nm %in% names(MASLD_DIR)) MASLD_DIR[[nm]] else "-"
  flag <- if (md != "-") (if (md == gd) "SAME" else "**DIFFERENT (disease-specific)**") else ""
  rows[[length(rows)+1]] <- data.frame(readout = nm, beta = round(b,3), p = round(p,4),
    GTEx_dir = gd, MASLD_dir = md, verdict = flag, stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(REALDIR, "gtex_healthy_sex.csv"), row.names = FALSE)
print(res, row.names = FALSE)
cat("\nKEY: for the MASLD headline readouts, does GTEx healthy liver show the SAME sex direction\n")
cat("(=> constitutional) or DIFFERENT/absent (=> MASLD-specific)?  Wrote gtex_healthy_sex.csv\n")
