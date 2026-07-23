## 12_liver_deconv.R ----------------------------------------------------------
## Liver-SPECIFIC deconvolution of the four MASLD bulk cohorts, using a reference
## built from a public single-cell liver atlas (Andrews 2024 by default). This
## complements 07_run_deconv.R (xCell / MCP-counter), which are trained on blood
## and tumour: a liver-matched reference is the natural next step and, unlike the
## bulk-reference tools, can carry a MAIT-cell cluster, giving an independent bulk
## estimate of the male-MAIT bias.
##
## Method: MuSiC (multi-subject single-cell reference -> bulk proportions).
## The reference is gated into liver cell types from canonical markers (the same
## MAIT / Treg / CD8 gates used across this study), then each cohort's symbol-level
## linear matrix (expr_sym_<GSE>.tsv from 06_prep_deconv.R) is deconvolved. The
## resulting per-sample fractions are then tested for a sex effect with the same
## fibrosis-adjusted mixed model used for the signature scores (03_maineffect.R).
##
## MAIT is a rare hepatic subset, so its deconvolved estimate is exploratory and
## is interpreted alongside the receptor-specific signature and the single-cell
## cohorts, not on its own.
## ===========================================================================
source("00_config.R")
suppressMessages({library(Seurat); library(Matrix)})

REF_RDS  <- Sys.getenv("MASLD_SC_REF",
              file.path(REALDIR, "single_cell/andrews2024/GSE243981_seurat.rds"))
REF_GROUP<- Sys.getenv("MASLD_SC_REF_GROUP", "healthy")   # restrict reference donors (NA/"" = all)
DONOR_COL<- "donor"
MAX_PER  <- 300     # cells kept per cell type per donor (reference size / speed)
set.seed(1)

## ---- liver reference cell types (hierarchical marker gates) ----------------
## T-lineage is resolved first (MAIT > Treg > CD8 > CD4), then non-T lineages,
## then non-immune parenchymal / stromal types, by strongest marker-set signal.
GATES_T <- list(
  MAIT = c("SLC4A10","TRAV1-2"),
  Treg = c("FOXP3","IL2RA","CTLA4"),
  CD8T = c("CD8A","CD8B"),
  CD4T = c("CD4","IL7R"))
GATES_NONT <- list(
  NK          = c("NCAM1","KLRF1","NKG7","GNLY"),
  Bcell       = c("CD19","MS4A1","CD79A"),
  Plasma      = c("MZB1","XBP1","DERL3"),
  cDC         = c("CD1C","CLEC9A","CLEC10A"),
  pDC         = c("LILRA4","IL3RA","GZMB"),
  MonoMac     = c("CD68","CD14","LYZ","CSF1R"),
  Neutrophil  = c("FCGR3B","CSF3R","S100A8"),
  Hepatocyte  = c("ALB","APOA1","TTR","APOB"),
  Endothelial = c("PECAM1","VWF","CLEC4G"),
  Cholangiocyte = c("EPCAM","KRT19","SOX9"),
  Stellate    = c("COL1A1","ACTA2","DCN","PDGFRB"))

## ---- load + subset the reference --------------------------------------------
message("Loading reference: ", REF_RDS)
if (!file.exists(REF_RDS)) stop("reference object not found; set MASLD_SC_REF")
ref <- readRDS(REF_RDS); DefaultAssay(ref) <- "RNA"
rmd <- ref@meta.data
if (nzchar(REF_GROUP) && REF_GROUP != "NA" && "group" %in% names(rmd))
  ref <- ref[, rmd$group == REF_GROUP]
rmd <- ref@meta.data
donor <- if (DONOR_COL %in% names(rmd)) as.character(rmd[[DONOR_COL]]) else rep("d1", ncol(ref))
cnt <- tryCatch(GetAssayData(ref, assay="RNA", layer="counts"),
                error=function(e) GetAssayData(ref, assay="RNA", slot="counts"))
genes <- rownames(cnt)
score <- function(gs) { gs <- gs[gs %in% genes]
  if (!length(gs)) return(rep(0, ncol(cnt)))
  Matrix::colMeans(cnt[gs, , drop=FALSE] > 0) }   # fraction of the marker set detected

## per-cell hierarchical assignment
isT <- score(c("CD3D","CD3E","CD3G")) > 0
Tsc  <- vapply(GATES_T,    score, numeric(ncol(cnt)))
NTsc <- vapply(GATES_NONT, score, numeric(ncol(cnt)))
lab <- character(ncol(cnt))
tpick  <- colnames(Tsc)[max.col(Tsc,  ties.method="first")]
ntpick <- colnames(NTsc)[max.col(NTsc, ties.method="first")]
lab[isT]  <- ifelse(apply(Tsc[isT,,drop=FALSE],1,max)  > 0, tpick[isT],  "CD4T")
lab[!isT] <- ifelse(apply(NTsc[!isT,,drop=FALSE],1,max) > 0, ntpick[!isT], NA)
keep <- !is.na(lab) & lab != ""
cat("Reference cells labelled:", sum(keep), "of", length(lab), "\n")
print(table(lab[keep]))

## cap cells per (type, donor) for a balanced, fast reference
idx <- which(keep)
grp <- paste(lab[idx], donor[idx])
sel <- unlist(lapply(split(idx, grp), function(ii) if (length(ii) > MAX_PER) sample(ii, MAX_PER) else ii))
refc <- cnt[, sel, drop=FALSE]; refl <- lab[sel]; refd <- donor[sel]
## drop ultra-rare types (< 20 cells) that MuSiC cannot estimate stably
tk <- names(table(refl))[table(refl) >= 20]
km <- refl %in% tk; refc <- refc[, km]; refl <- refl[km]; refd <- refd[km]
cat("Reference used:", ncol(refc), "cells,", length(unique(refl)), "types,",
    length(unique(refd)), "donors\n")

## ---- MuSiC deconvolution per cohort -----------------------------------------
suppressMessages({library(MuSiC); library(SingleCellExperiment)})
sce <- SingleCellExperiment(assays = list(counts = as.matrix(refc)),
                            colData = DataFrame(celltype = refl, donor = refd))
props <- list()
for (g in COHORT_IDS) {
  fp <- file.path(REALDIR, paste0("expr_sym_", g, ".tsv"))
  if (!file.exists(fp)) { cat("missing bulk matrix:", fp, "\n"); next }
  B <- as.data.frame(data.table::fread(fp)); rn <- as.character(B[[1]]); B[[1]] <- NULL
  B <- as.matrix(B); rownames(B) <- rn; storage.mode(B) <- "double"
  est <- tryCatch(
    MuSiC::music_prop(bulk.mtx = B, sc.sce = sce, clusters = "celltype",
                      samples = "donor", verbose = FALSE),
    error = function(e) { cat("MuSiC error", g, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(est)) next
  P <- as.data.frame(est$Est.prop.weighted)
  P$sample <- rownames(P); P$cohort <- g
  props[[g]] <- P
  cat(g, "deconvolved:", nrow(P), "samples\n")
}
if (!length(props)) stop("no cohorts deconvolved")
allP <- do.call(function(...) rbind(..., make.row.names=FALSE),
                lapply(props, function(d) { ct <- setdiff(names(d), c("sample","cohort"))
                  d[, c("sample","cohort", ct)] }))
write.csv(allP, file.path(REALDIR, "liver_deconv_props.csv"), row.names=FALSE)
cat("Wrote liver_deconv_props.csv\n")

## ---- sex effect on the deconvolved fractions --------------------------------
suppressMessages({library(lme4); library(lmerTest)})
mtx <- read.csv(file.path(REALDIR, "analysis_matrix.csv"), stringsAsFactors=FALSE)
key <- intersect(c("sample","gsm","geo_accession"), names(mtx))[1]
if (is.na(key)) key <- names(mtx)[1]
mtx$fib_ord <- code_fibrosis(mtx$fibrosis_stage)
mtx$sexM    <- as.numeric(mtx$sex_assigned == "M")
M <- merge(allP, mtx[, c(key, "cohort", "sex_assigned", "sexM", "fib_ord")],
           by.x = c("sample","cohort"), by.y = c(key, "cohort"))
M <- M[M$sex_assigned %in% c("F","M") & !is.na(M$fib_ord), , drop=FALSE]
celltypes <- setdiff(names(allP), c("sample","cohort"))

## ---- resolution diagnostics -----------------------------------------------
## Whole-liver bulk is ~90% hepatocyte, so rare subsets (MAIT, Treg, CD8, CD4,
## B cells) are pushed toward zero and are NOT reliably deconvolvable. Report
## which types carry real signal, and only model those; the rest are flagged
## "below whole-tissue bulk resolution" rather than returned as spurious NA/0.
diag <- data.frame(celltype = celltypes,
  mean_prop = round(colMeans(allP[, celltypes], na.rm=TRUE), 4),
  pct_nearzero = round(100*colMeans(allP[, celltypes] < 1e-4, na.rm=TRUE), 1),
  sd_prop = round(apply(allP[, celltypes], 2, sd, na.rm=TRUE), 4), row.names=NULL)
diag <- diag[order(-diag$mean_prop), ]
cat("\n== deconvolved proportions: mean / % near-zero / sd (per cell type) ==\n")
print(diag, row.names=FALSE)
resolvable <- diag$celltype[diag$pct_nearzero < 80 & diag$sd_prop > 1e-4]
below      <- setdiff(celltypes, resolvable)
cat("\nModelled (resolvable):", paste(resolvable, collapse=", "), "\n")
if (length(below)) cat("Below whole-liver bulk resolution (not modelled):",
                        paste(below, collapse=", "), "\n")

## within-cohort z-score each resolvable fraction (comparable betas across studies)
for (ct in resolvable) M[[ct]] <- ave(M[[ct]], M$cohort, FUN=function(x) zscore(x))

## mixed model with a fixed-effect (cohort as covariate) fallback if it is singular
fit1 <- function(ct) {
  d <- M[stats::complete.cases(M[, c(ct,"fib_ord","sexM","cohort")]), ]
  m <- tryCatch(lmerTest::lmer(as.formula(sprintf("`%s` ~ sexM + fib_ord + (1|cohort)", ct)),
        d, REML=TRUE, control=lmerControl(calc.derivs=FALSE)), error=function(e) NULL)
  co <- if (!is.null(m)) tryCatch(summary(m)$coefficients, error=function(e) NULL) else NULL
  if (!is.null(co) && "sexM" %in% rownames(co) && !is.na(co["sexM","Std. Error"]))
    return(c(beta=co["sexM","Estimate"], se=co["sexM","Std. Error"], p=co["sexM","Pr(>|t|)"], n=nrow(d), model=1))
  fit <- lm(as.formula(sprintf("`%s` ~ sexM + fib_ord + cohort", ct)), d)  # fixed-effect fallback
  c2 <- summary(fit)$coefficients
  c(beta=c2["sexM","Estimate"], se=c2["sexM","Std. Error"], p=c2["sexM","Pr(>|t|)"], n=nrow(d), model=2)
}
res <- if (length(resolvable))
  do.call(rbind, lapply(resolvable, function(ct) data.frame(celltype=ct, t(fit1(ct))))) else NULL
if (!is.null(res)) {
  res$fdr <- p.adjust(res$p, method="BH")
  res$direction <- ifelse(res$beta>0, "male", "female")
  res$model <- ifelse(res$model==1, "mixed", "fixed-fallback")
  res <- res[order(res$fdr), ]
}
## append the below-resolution types for transparency
if (length(below)) {
  br <- data.frame(celltype=below, beta=NA, se=NA, p=NA, n=NA, model="below_resolution",
                   fdr=NA, direction=NA)
  res <- if (is.null(res)) br else rbind(res[, names(br)], br)
}
write.csv(res, file.path(REALDIR, "liver_deconv_sex.csv"), row.names=FALSE)
cat("\n==== sex effect on liver-reference deconvolved fractions (beta>0 = male) ====\n")
print(res, row.names=FALSE, digits=3)
cat("\nNote: rare hepatic T-cell subsets (incl. MAIT) are typically below the",
    "\nresolution of whole-tissue bulk deconvolution against a hepatocyte-dominated",
    "\nreference; this is expected and is why MAIT rests on the receptor-specific",
    "\nsignature and single-cell cohorts rather than bulk deconvolution.\n")
cat("Wrote liver_deconv_sex.csv and liver_deconv_props.csv to", REALDIR, "\n")
