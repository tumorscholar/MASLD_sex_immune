## sc_validate_andrews.R ------------------------------------------------------
## Single-cell validation of the sex x hepatic-immune findings in Andrews 2024
## (GSE243981). Tests the four headline directions at single-cell resolution:
##   MAIT higher in MALES ;  Treg / CD8 / cDC higher in FEMALES
## Reports OVERALL and WITHIN each disease group (the healthy donors are the
## clean test of the CONSTITUTIONAL male-MAIT prediction).
##
## Method: pseudobulk each donor -> assign sex (XIST vs Y panel); gate the four
## populations per cell from markers; per-donor fractions -> Wilcoxon by sex.
## ===========================================================================
OBJ_PATH   <- file.path(Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta"), "single_cell/andrews2024/GSE243981_seurat.rds")
DONOR_COL  <- "donor"    # per the build script
GROUP_COL  <- "group"    # healthy / PSC / PBC  (set NA if none)
SEX_COL    <- NA         # Andrews doesn't record sex in GEO -> assign from expression
OUTDIR     <- file.path(Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta"), "single_cell/andrews2024/out")
## ===========================================================================
DATASET <- "Andrews2024"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
suppressMessages({library(Seurat); library(Matrix); library(ggplot2)})

FEMALE <- "XIST"
MALE   <- c("RPS4Y1","DDX3Y","EIF1AY","UTY","KDM5D","USP9Y","NLGN4Y","ZFY","TXLNGY")
MARK <- list(
 Tcell = c("CD3D","CD3E","CD3G"),
 MAIT  = c("SLC4A10","TRAV1-2"),
 CD8   = c("CD8A","CD8B"),
 Treg  = c("FOXP3"),
 cDC   = c("CD1C","CLEC9A","CLEC10A","FCER1A"),
 pDC   = c("LILRA4","IL3RA","GZMB"),
 immune= c("PTPRC"))

message("Loading ", OBJ_PATH)
obj <- readRDS(OBJ_PATH); DefaultAssay(obj) <- "RNA"
md  <- obj@meta.data
stopifnot(DONOR_COL %in% colnames(md))
donor  <- as.character(md[[DONOR_COL]]); donors <- sort(unique(donor))
group  <- if (!is.na(GROUP_COL) && GROUP_COL %in% colnames(md)) as.character(md[[GROUP_COL]]) else rep("all", nrow(md))
donor_group <- tapply(group, donor, function(x) x[1])[donors]
cat("Donors:", length(donors), "| cells:", ncol(obj), "\n")

cnt <- tryCatch(GetAssayData(obj, assay="RNA", layer="counts"),
                error=function(e) GetAssayData(obj, assay="RNA", slot="counts"))
genes <- rownames(cnt); has <- function(g) g %in% genes
gvec <- function(g) if (has(g)) as.numeric(cnt[g, ]) else rep(0, ncol(cnt))
anyMark <- function(gs){ gs<-gs[vapply(gs,has,logical(1))]; if(!length(gs)) return(rep(FALSE,ncol(cnt)))
Reduce(`|`, lapply(gs, function(g) gvec(g)>0)) }
for (nm in names(MARK)) { miss <- MARK[[nm]][!vapply(MARK[[nm]],has,logical(1))]
if (length(miss)) cat("  note: missing", nm, "markers:", paste(miss,collapse=","), "\n") }

## 1. per-donor pseudobulk sex assignment
zscore <- function(x){ m<-mean(x,na.rm=TRUE); s<-sqrt(mean((x-m)^2,na.rm=TRUE)); if(is.na(s)||s==0) x*0 else (x-m)/s }
sexgenes <- intersect(c(FEMALE,MALE), genes)
pb <- t(vapply(donors, function(d) Matrix::rowSums(cnt[sexgenes, donor==d, drop=FALSE]), numeric(length(sexgenes))))
colnames(pb) <- sexgenes
libsize <- vapply(donors, function(d) sum(cnt[, donor==d]), numeric(1))
pb <- log2(sweep(pb,1,libsize,"/")*1e6 + 1)
zX <- if (FEMALE %in% colnames(pb)) zscore(pb[,FEMALE]) else rep(0,length(donors))
yy <- intersect(MALE, colnames(pb))
zY <- if (length(yy)) zscore(rowMeans(pb[,yy,drop=FALSE])) else rep(0,length(donors))
sex_assigned <- ifelse(zY-zX > 0.5, "M", ifelse(zY-zX < -0.5, "F", "Ambiguous")); names(sex_assigned) <- donors
cat("\nAssigned sex:", sum(sex_assigned=="M"),"M /", sum(sex_assigned=="F"),"F /",
    sum(sex_assigned=="Ambiguous"),"ambiguous\n")
if (!is.na(SEX_COL) && SEX_COL %in% colnames(md)) {
 rec <- tapply(as.character(md[[SEX_COL]]), donor, function(x) tolower(x[1]))[donors]
 rec <- ifelse(startsWith(rec,"m"),"M", ifelse(startsWith(rec,"f"),"F",NA))
 ok <- !is.na(rec) & sex_assigned!="Ambiguous"
 if (any(ok)) cat(sprintf("Concordance with recorded sex: %.1f%% (%d donors)\n",
                          100*mean(rec[ok]==sex_assigned[ok]), sum(ok)))
}

## 2. per-cell gates
is_T   <- anyMark(MARK$Tcell)
is_MAIT<- is_T & anyMark(MARK$MAIT)
is_CD8 <- is_T & anyMark(MARK$CD8)
is_Treg<- is_T & anyMark(MARK$Treg)
is_imm <- anyMark(MARK$immune) | is_T
is_cDC <- anyMark(MARK$cDC) & !anyMark(MARK$pDC)

## 3. per-donor fractions
frac <- function(numer,denom){ n<-tapply(numer,donor,sum)[donors]; d<-tapply(denom,donor,sum)[donors]; ifelse(d>0,100*n/d,NA) }
df <- data.frame(donor=donors, group=donor_group[donors], sex=sex_assigned[donors],
                 n_cells=as.integer(table(donor)[donors]), n_Tcells=as.integer(tapply(is_T,donor,sum)[donors]),
                 n_MAIT=as.integer(tapply(is_MAIT,donor,sum)[donors]),
                 MAIT_pctT=frac(is_MAIT,is_T), Treg_pctT=frac(is_Treg,is_T),
                 CD8_pctT=frac(is_CD8,is_T), cDC_pctImm=frac(is_cDC,is_imm), row.names=NULL)
write.csv(df, file.path(OUTDIR, paste0(DATASET,"_per_donor_fractions.csv")), row.names=FALSE)
cat("\nMAIT cells per donor (median):", median(df$n_MAIT,na.rm=TRUE),
    "| donors with >=10 MAIT:", sum(df$n_MAIT>=10,na.rm=TRUE), "/", nrow(df), "\n")

## 4. test overall + within each group
EXPECT <- c(MAIT_pctT="M", Treg_pctT="F", CD8_pctT="F", cDC_pctImm="F")
report <- function(d2, label) {
 d2 <- d2[d2$sex %in% c("M","F"), ]
 if (sum(d2$sex=="M")<2 || sum(d2$sex=="F")<2) { cat("\n[",label,"] too few per sex (M=",sum(d2$sex=="M"),
                                                     " F=",sum(d2$sex=="F"),") - skipped\n",sep=""); return(NULL) }
 res <- do.call(rbind, lapply(names(EXPECT), function(v){
  m<-d2[[v]][d2$sex=="M"]; f<-d2[[v]][d2$sex=="F"]
  p<-tryCatch(wilcox.test(m,f)$p.value, error=function(e) NA)
  hi<-if(median(m,na.rm=TRUE)>median(f,na.rm=TRUE)) "M" else "F"
  data.frame(readout=v, median_M=round(median(m,na.rm=TRUE),2), median_F=round(median(f,na.rm=TRUE),2),
             higher_in=hi, predicted=EXPECT[[v]], matches=ifelse(hi==EXPECT[[v]],"YES","no"), p=signif(p,3), row.names=NULL)}))
 cat("\n==== ",label," (M=",sum(d2$sex=="M"),", F=",sum(d2$sex=="F")," donors) ====\n",sep="")
 print(res, row.names=FALSE); cat("Directions reproduced:", sum(res$matches=="YES"),"/ 4\n")
 res$subset<-label; res }
allres <- report(df, "ALL donors")
for (g in setdiff(unique(na.omit(df$group)), "all")) allres <- rbind(allres, report(df[df$group==g,], g))
if (!is.null(allres)) write.csv(allres, file.path(OUTDIR, paste0(DATASET,"_sex_validation.csv")), row.names=FALSE)

## plot (healthy donors, or all if no groups)
plotset <- if ("healthy" %in% df$group) df[df$group=="healthy" & df$sex %in% c("M","F"),] else df[df$sex %in% c("M","F"),]
long <- do.call(rbind, lapply(names(EXPECT), function(v) data.frame(readout=v, sex=plotset$sex, value=plotset[[v]])))
long$readout <- factor(long$readout, levels=names(EXPECT), labels=c("MAIT (%T)","Treg (%T)","CD8 (%T)","cDC (%imm)"))
p <- ggplot(long, aes(sex,value,fill=sex)) + geom_boxplot(outlier.shape=NA,alpha=.5) +
 geom_jitter(width=.15,size=1.6,alpha=.8) + facet_wrap(~readout,scales="free_y",nrow=1) +
 scale_fill_manual(values=c(M="#2c7fb8",F="#d95f0e")) +
 labs(x=NULL,y="per-donor fraction (%)", title=paste0(DATASET,": hepatic immune fractions by sex")) +
 theme_classic(base_size=10) + theme(legend.position="none")
ggsave(file.path(OUTDIR, paste0(DATASET,"_sex_validation.png")), p, width=9, height=3, dpi=300)
cat("\nWrote fractions, results and plot to", OUTDIR, "\n")
