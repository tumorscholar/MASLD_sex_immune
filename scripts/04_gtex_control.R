## 04_gtex_control.R ----------------------------------------------------------
## Disease-specificity control. Score the identical immune readouts in GTEx
## disease-free liver and test the sex effect, to ask whether each MASLD sex bias
## is CONSTITUTIONAL (present in healthy liver too) or DISEASE-EMERGENT.
##
## The sex effect is estimated under three increasingly strict models, so the
## conclusion does not rest on a minimally-adjusted comparison:
##   A) age-only            readout ~ sexM + age
##   B) covariate-adjusted  readout ~ sexM + age + RIN + ischaemic time + centre
##   C) sex-balanced        model B on K balanced male/female subsamples (GTEx
##      liver is male-skewed); reports median beta + direction consistency.
## RIN = SMRIN, ischaemic time = SMTSISCH, centre = SMCENTER (SampleAttributes).
## ---------------------------------------------------------------------------
source("00_config.R")
suppressMessages({library(data.table); library(singscore)})
set.seed(1)

## Score a signature with the Bioconductor singscore package (identical to
## 02_build_matrix.R). Filters gene IDs to those present in the ranked matrix and
## returns the per-sample TotalScore; NULL if the up-set has no genes present.
score_sig <- function(ranks, up_ids, dn_ids = NULL) {
  up_ids <- unique(up_ids[!is.na(up_ids) & up_ids %in% rownames(ranks)])
  if (!length(up_ids)) return(NULL)
  dn_ids <- if (!is.null(dn_ids)) unique(dn_ids[!is.na(dn_ids) & dn_ids %in% rownames(ranks)]) else NULL
  sc <- if (!is.null(dn_ids) && length(dn_ids))
    singscore::simpleScore(ranks, upSet = up_ids, downSet = dn_ids)
  else
    singscore::simpleScore(ranks, upSet = up_ids)
  sc$TotalScore
}
K_BAL <- 100   # balanced subsamples for model C

BASE <- "https://storage.googleapis.com/adult-gtex"
FILES <- c(
  attr = paste0(BASE, "/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
  phen = paste0(BASE, "/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"),
  tpm  = paste0(BASE, "/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
)
fetch <- function(k) {
  p <- file.path(SCRATCH, basename(FILES[k]))
  if (!file.exists(p)) { message("  downloading ", k); download.file(FILES[k], p, mode="wb", quiet=TRUE) }
  p
}

GT_CELLTYPE <- CELLTYPE
GT_CELLTYPE$ct_MAITspec    <- c("SLC4A10","TRAV1-2")
GT_CELLTYPE$ct_MAITpromisc <- c("KLRB1","RORC","ZBTB16")
GT_STATE <- STATE[c("st_CD8_cytotox","st_Th1","st_Th17","st_cytotoxCD4","st_senescence")]
MASLD_DIR <- c(ct_MAIT="MALE", ct_MAITspec="MALE", ct_Treg="female",
               ct_CD8T="female", ct_DC="female", st_Th1="female")

## ---- annotations (keep the technical covariates for liver samples) ----
message("Fetching GTEx annotations...")
attr <- fread(fetch("attr"), sep = "\t")
phen <- fread(fetch("phen"), sep = "\t")
liverA <- attr[SMTSD == "Liver", .(SAMPID, SMRIN, SMTSISCH, SMCENTER)]
liver  <- liverA$SAMPID
cat("GTEx liver samples:", length(liver), "\n")
setkey(phen, SUBJID)
sid <- function(s) vapply(strsplit(s, "-"), function(x) paste(x[1:2], collapse="-"), character(1))

## ---- stream the gene-TPM gct, liver columns only (memory-safe) ----
message("Reading GTEx gene TPM (streaming liver columns only)...")
con <- gzfile(fetch("tpm"), "rt")
readLines(con, 2); header <- strsplit(readLines(con, 1), "\t")[[1]]
di <- which(header == "Description"); liv <- liver[liver %in% header]; li <- match(liv, header)
genes <- character(0); mat <- vector("list", 0); k <- 0
repeat {
  chunk <- readLines(con, 2000); if (length(chunk) == 0) break
  sp <- strsplit(chunk, "\t")
  genes <- c(genes, vapply(sp, function(p) p[di], character(1)))
  for (p in sp) { k <- k + 1; mat[[k]] <- as.numeric(p[li]) }
}
close(con)
M <- do.call(rbind, mat); rownames(M) <- genes; colnames(M) <- liv
M <- M[!duplicated(rownames(M)), , drop = FALSE]; M <- log2(pmax(M, 0) + 1)
cat("Expression matrix:", nrow(M), "x", ncol(M), "\n")

## ---- score ----
ranks <- singscore::rankGenes(M, tiesMethod = "average"); score <- list()
for (nm in names(GT_CELLTYPE)) { s <- score_sig(ranks, GT_CELLTYPE[[nm]]); if(!is.null(s)) score[[nm]] <- s }
for (nm in names(GT_STATE))    { s <- score_sig(ranks, GT_STATE[[nm]]$up, GT_STATE[[nm]]$down); if(!is.null(s)) score[[nm]] <- s }
## singscore::simpleScore()$TotalScore returns an unnamed vector, so restore the
## GTEx sample IDs as rownames (they are in colnames(M) order) before deriving
## SAMPID – otherwise the merge with donor/technical metadata finds no matches.
S <- as.data.frame(score); rownames(S) <- colnames(M); S$SAMPID <- rownames(S)

## ---- donor + technical metadata ----
S$subj <- sid(S$SAMPID); S$SEX <- phen[S$subj, SEX]; S$AGE <- phen[S$subj, AGE]
S$sexM <- as.numeric(S$SEX == 1)
agemid <- function(a) vapply(a, function(v){ z<-strsplit(v,"-")[[1]]; if(length(z)==2)(as.numeric(z[1])+as.numeric(z[2]))/2 else NA_real_ }, numeric(1))
S$age <- agemid(S$AGE)
S <- merge(S, as.data.frame(liverA), by = "SAMPID", all.x = TRUE)
S$SMRIN <- as.numeric(S$SMRIN); S$SMTSISCH <- as.numeric(S$SMTSISCH); S$SMCENTER <- factor(S$SMCENTER)
S <- S[!is.na(S$sexM), , drop = FALSE]
cat(sprintf("GTEx liver donors: M=%d F=%d | RIN available=%d, ischaemic-time=%d\n",
            sum(S$sexM==1), sum(S$sexM==0), sum(!is.na(S$SMRIN)), sum(!is.na(S$SMTSISCH))))

## ---- sex-effect estimator under an arbitrary right-hand side ----
fit_beta <- function(d, nm, rhs) {
  vars <- c(nm, all.vars(as.formula(paste("~", rhs))))
  d <- d[stats::complete.cases(d[, intersect(vars, names(d))]), , drop = FALSE]
  if (nrow(d) < 10) return(c(beta=NA, p=NA, n=nrow(d)))
  d[[nm]] <- zscore(d[[nm]])
  f <- tryCatch(lm(as.formula(paste0(nm, " ~ ", rhs)), d), error = function(e) NULL)
  if (is.null(f)) return(c(beta=NA, p=NA, n=nrow(d)))
  co <- summary(f)$coefficients
  if (!"sexM" %in% rownames(co)) return(c(beta=NA, p=NA, n=nrow(d)))
  c(beta=co["sexM","Estimate"], p=co["sexM","Pr(>|t|)"], n=nrow(d))
}
COV_RHS <- "sexM + age + SMRIN + SMTSISCH + SMCENTER"
balanced <- function(d, nm, K = K_BAL) {
  vars <- c(nm, "sexM","age","SMRIN","SMTSISCH","SMCENTER")
  d <- d[stats::complete.cases(d[, vars]), , drop = FALSE]
  Fi <- which(d$sexM==0); Mi <- which(d$sexM==1); nF <- length(Fi)
  if (nF < 5 || length(Mi) < 5) return(c(beta=NA, consist=NA))
  bs <- numeric(0)
  for (i in 1:K) { samp <- d[c(Fi, sample(Mi, nF)), , drop=FALSE]
  b <- fit_beta(samp, nm, COV_RHS)["beta"]; if (!is.na(b)) bs <- c(bs, b) }
  md <- median(bs, na.rm=TRUE); c(beta=md, consist=mean(sign(bs)==sign(md)))
}

## ---- run all three models per readout ----
cat("\n==== SEX effect in DISEASE-FREE GTEx liver – three models ====\n")
readout_names <- c(names(GT_CELLTYPE), names(GT_STATE)); rows <- list()
for (nm in readout_names) {
  if (!nm %in% names(S)) next
  a <- fit_beta(S, nm, "sexM + age")                 # A age-only
  b <- fit_beta(S, nm, COV_RHS)                       # B covariate-adjusted
  cc <- balanced(S, nm)                               # C sex-balanced (median of K)
  dir <- function(x) if (is.na(x)) "-" else if (x>0) "MALE" else "female"
  md  <- if (nm %in% names(MASLD_DIR)) MASLD_DIR[[nm]] else "-"
  gdc <- dir(b["beta"])                               # covariate-adjusted direction = primary
  verdict <- if (md=="-") "" else if (gdc==md) "SAME -> constitutional" else "DIFFERENT -> disease-specific"
  rows[[length(rows)+1]] <- data.frame(readout=nm, MASLD_dir=md,
                                       beta=round(b["beta"],3),                       # primary = covariate-adjusted (used by Fig 5)
                                       beta_ageonly=round(a["beta"],3), dir_ageonly=dir(a["beta"]),
                                       beta_cov=round(b["beta"],3), p_cov=round(b["p"],4), dir_cov=gdc, n_cov=b["n"],
                                       beta_bal=round(cc["beta"],3), dir_bal=dir(cc["beta"]), bal_consistency=round(cc["consist"],2),
                                       verdict=verdict, row.names=NULL, stringsAsFactors=FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(REALDIR, "gtex_healthy_sex.csv"), row.names = FALSE)
print(res, row.names = FALSE)
cat("\nColumns: beta>0 = higher in MALES. dir_cov (covariate-adjusted) is the primary read;",
    "\ndir_bal + bal_consistency show robustness to the male sex skew (fraction of", K_BAL,
    "\nbalanced subsamples agreeing with the median direction). Wrote gtex_healthy_sex.csv.\n")
cat("\nKEY: for the MASLD headline readouts, is the sex direction preserved under the FULL\n")
cat("covariate + sex-balanced model? MAIT staying MALE = robust constitutional bias; Treg/CD8/DC\n")
cat("staying female-absent/weak = robust disease-emergence. This addresses the age-only caveat.\n")