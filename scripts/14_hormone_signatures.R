## 14_hormone_signatures.R ----------------------------------------------------
## Mechanistic support for the sex-hormone hypothesis in the Discussion. Scores
## androgen- and oestrogen-response transcriptional programmes per sample and asks
## two questions in the bulk cohorts:
##   (1) do these hormone-response signatures themselves differ by sex?
##   (2) do they track the sex-divergent immune axes, i.e. does a higher
##       androgen-response programme go with a lower female-biased immune
##       programme (Treg / CD8 / cDC) and does oestrogen-response relate to MAIT?
## This turns the hormone hypothesis (currently cited, not shown) into a
## computational observation, without any new data.
##
## Signatures: MSigDB Hallmark ANDROGEN_RESPONSE, ESTROGEN_RESPONSE_EARLY/LATE.
## Scoring: the same rank-based single-sample method (singscore) used throughout.
## ===========================================================================
source("00_config.R")
suppressMessages({library(data.table); library(lme4); library(lmerTest)})
if (!requireNamespace("msigdbr", quietly=TRUE))
  stop("msigdbr needed: run install_packages.R")
suppressMessages(library(msigdbr))

## ---- hormone-response gene sets --------------------------------------------
H <- msigdbr(species="Homo sapiens", category="H")
gs <- function(nm) unique(H$gene_symbol[H$gs_name == nm])
SETS <- list(
  Androgen_response = gs("HALLMARK_ANDROGEN_RESPONSE"),
  Estrogen_early    = gs("HALLMARK_ESTROGEN_RESPONSE_EARLY"),
  Estrogen_late     = gs("HALLMARK_ESTROGEN_RESPONSE_LATE"))
for (nm in names(SETS)) cat(sprintf("  %-18s %d genes\n", nm, length(SETS[[nm]])))

## ---- score per sample in each cohort (rank-based singscore) -----------------
scores <- list()
for (g in COHORT_IDS) {
  fp <- file.path(REALDIR, paste0("expr_sym_", g, ".tsv"))
  if (!file.exists(fp)) { cat("missing", fp, "\n"); next }
  M <- as.data.frame(fread(fp)); rn <- as.character(M[[1]]); M[[1]] <- NULL
  M <- as.matrix(M); rownames(M) <- rn; storage.mode(M) <- "double"
  M <- M[rowSums(!is.na(M)) > 0, , drop=FALSE]
  ng <- nrow(M); ranks <- rank_matrix(M)                 # within-sample ranks (scale-invariant)
  s <- sapply(names(SETS), function(nm) {
    v <- singscore_custom(ranks, ng, SETS[[nm]]); if (is.null(v)) rep(NA, ncol(M)) else v })
  d <- as.data.frame(s); d$sample <- colnames(M); d$cohort <- g
  scores[[g]] <- d
  cat(g, "scored:", ncol(M), "samples\n")
}
S <- as.data.frame(rbindlist(scores, fill=TRUE))
write.csv(S, file.path(REALDIR, "hormone_scores.csv"), row.names=FALSE)

## ---- merge with the analysis matrix (sex, fibrosis, immune readouts) --------
mtx <- read.csv(file.path(REALDIR, "analysis_matrix.csv"), stringsAsFactors=FALSE)
key <- intersect(c("sample","gsm","geo_accession"), names(mtx))[1]
if (is.na(key)) key <- names(mtx)[1]
mtx$fib_ord <- code_fibrosis(mtx$fibrosis_stage)
mtx$sexM    <- as.numeric(mtx$sex_assigned == "M")
D <- merge(S, mtx, by.x=c("sample","cohort"), by.y=c(key,"cohort"))
D <- D[D$sex_assigned %in% c("F","M") & !is.na(D$fib_ord), , drop=FALSE]
HORM <- names(SETS)
## within-cohort z-score the hormone scores (comparable across cohorts)
for (h in HORM) D[[h]] <- ave(D[[h]], D$cohort, FUN=function(x) zscore(x))

## ---- (1) sex effect on the hormone-response signatures ----------------------
sex_on <- function(h) {
  d <- D[stats::complete.cases(D[, c(h,"fib_ord","sexM","cohort")]), ]
  out <- tryCatch({ m <- lmerTest::lmer(as.formula(sprintf("`%s` ~ sexM + fib_ord + (1|cohort)", h)),
      d, REML=TRUE, control=lmerControl(calc.derivs=FALSE))
    co <- summary(m)$coefficients; c(beta=co["sexM","Estimate"], se=co["sexM","Std. Error"], p=co["sexM","Pr(>|t|)"]) },
    error=function(e) c(beta=NA, se=NA, p=NA))
  data.frame(signature=h, beta_sexM=out["beta"], p=out["p"],
             direction=ifelse(is.na(out["beta"]),NA,ifelse(out["beta"]>0,"male","female")), row.names=NULL)
}
t_sex <- do.call(rbind, lapply(HORM, sex_on))
t_sex$fdr <- p.adjust(t_sex$p, method="BH")
write.csv(t_sex, file.path(REALDIR, "hormone_sex_effect.csv"), row.names=FALSE)
cat("\n==== (1) sex effect on hormone-response signatures (beta>0 = higher in males) ====\n")
print(t_sex, row.names=FALSE, digits=3)

## ---- (2) do the hormone programmes track the immune axes? -------------------
## For each headline immune readout, does the hormone signature predict it beyond
## sex and fibrosis? (coefficient on the hormone score in readout ~ horm + sexM + fib)
READS <- intersect(c("ct_MAIT","ct_MAITspec","ct_Treg","ct_CD8T","ct_DC","st_Th1"), names(D))
for (r in READS) D[[r]] <- ave(D[[r]], D$cohort, FUN=function(x) zscore(x))
assoc <- list()
for (r in READS) for (h in HORM) {
  d <- D[stats::complete.cases(D[, c(r,h,"fib_ord","sexM","cohort")]), ]
  out <- tryCatch({ m <- lmerTest::lmer(as.formula(sprintf("`%s` ~ `%s` + sexM + fib_ord + (1|cohort)", r, h)),
      d, REML=TRUE, control=lmerControl(calc.derivs=FALSE))
    co <- summary(m)$coefficients; c(b=co[h,"Estimate"], p=co[h,"Pr(>|t|)"]) },
    error=function(e) c(b=NA, p=NA))
  assoc[[length(assoc)+1]] <- data.frame(readout=r, hormone=h, beta=out["b"], p=out["p"], row.names=NULL)
}
A <- do.call(rbind, assoc); A$fdr <- p.adjust(A$p, method="BH")
A <- A[order(A$fdr), ]
write.csv(A, file.path(REALDIR, "hormone_immune_assoc.csv"), row.names=FALSE)
cat("\n==== (2) hormone-response vs immune readouts (partialling out sex + fibrosis) ====\n")
cat("beta>0 = readout rises with the hormone programme; beta<0 = falls.\n")
print(A, row.names=FALSE, digits=3)
cat("\nExpected pattern if the hypothesis holds: androgen-response tracks NEGATIVELY",
    "\nwith the female-biased programme (Treg/CD8/cDC); oestrogen-response relates to MAIT.\n")
cat("Wrote hormone_scores.csv, hormone_sex_effect.csv, hormone_immune_assoc.csv to", REALDIR, "\n")
