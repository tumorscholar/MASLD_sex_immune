## sc_mait_functional.R -------------------------------------------------------
## Functional comparison of MALE vs FEMALE hepatic MAIT cells at the same stage
## (Andrews healthy liver). Asks: beyond being more frequent in men, are male MAIT
## cells transcriptionally different (e.g. more cytotoxic / pro-inflammatory)?
##
## Design:
##   1. gate MAIT (T cell & SLC4A10/TRAV1-2+) in healthy donors
##   2. PSEUDOBULK MAIT per donor (sum counts) -> no per-cell pseudoreplication
##   3. male-vs-female DE with edgeR quasi-likelihood
##   4. fgsea on the ranked genes (Reactome + Hallmark, via msigdbr)
##   5. functional-state scores per donor (cytotoxicity / cytokine / exhaustion /
##      activation) compared by sex
## Outputs DE + fgsea + state tables and figures (PDF+TIFF+PNG).
## ===========================================================================
OBJ_PATH  <- file.path(Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta"), "single_cell/andrews2024/GSE243981_seurat.rds")
GROUP_KEEP<- "healthy"     # same-stage comparison (constitutional MAIT)
MIN_MAIT  <- 10            # donors need >= this many MAIT cells to contribute
OUTDIR    <- file.path(Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta"), "single_cell/andrews2024/out")
FIGDIR    <- file.path(Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta"), "figures")   # Supp Fig S3 goes here, with the other figures
## ===========================================================================
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
dir.create(FIGDIR, showWarnings=FALSE, recursive=TRUE)
need <- function(p, bioc=FALSE) if (!requireNamespace(p, quietly=TRUE)) {
  if (bioc) BiocManager::install(p, update=FALSE, ask=FALSE) else install.packages(p) }
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (p in c("edgeR","limma","fgsea")) need(p, bioc=TRUE)
for (p in c("msigdbr","ggplot2")) need(p)
suppressMessages({library(Seurat); library(Matrix); library(edgeR); library(fgsea); library(msigdbr); library(ggplot2)})

save_fig <- function(draw, name, w, h) {
  render <- function() if (inherits(draw,"ggplot")) print(draw) else draw()
  pdf (file.path(FIGDIR, paste0(name,".pdf")),  width=w, height=h);                                       render(); dev.off()
  tiff(file.path(FIGDIR, paste0(name,".tiff")), width=w, height=h, units="in", res=600, compression="lzw"); render(); dev.off()
  cat("wrote", name, "-> ", FIGDIR, "\n")
}
MALE <- "#2c7fb8"; FEM <- "#d95f0e"
zscore <- function(x){ m<-mean(x,na.rm=TRUE); s<-sd(x,na.rm=TRUE); if(is.na(s)||s==0) x*0 else (x-m)/s }

## ---- load + restrict to same-stage donors ----
obj <- readRDS(OBJ_PATH); DefaultAssay(obj) <- "RNA"; md <- obj@meta.data
if (!is.na(GROUP_KEEP) && "group" %in% names(md)) obj <- obj[, md$group == GROUP_KEEP]
md <- obj@meta.data; donor <- as.character(md$donor)
cnt <- tryCatch(GetAssayData(obj, assay="RNA", layer="counts"), error=function(e) GetAssayData(obj, assay="RNA", slot="counts"))
genes <- rownames(cnt); has <- function(g) g %in% genes
gv <- function(g) if (has(g)) as.numeric(cnt[g,]) else rep(0, ncol(cnt))
anyMark <- function(gs){ gs<-gs[vapply(gs,has,logical(1))]; if(!length(gs)) return(rep(FALSE,ncol(cnt))); Reduce(`|`, lapply(gs,function(g) gv(g)>0)) }

## ---- gate MAIT ----
is_T    <- anyMark(c("CD3D","CD3E","CD3G"))
is_MAIT <- is_T & anyMark(c("SLC4A10","TRAV1-2"))
cat("MAIT cells:", sum(is_MAIT), "in", ncol(obj), "cells (", GROUP_KEEP, ")\n")

## ---- sex per donor (pseudobulk XIST vs Y on all cells) ----
FEMALE<-"XIST"; MALEG<-c("RPS4Y1","DDX3Y","EIF1AY","UTY","KDM5D","USP9Y","NLGN4Y","ZFY","TXLNGY")
pts <- sort(unique(donor)); sg <- intersect(c(FEMALE,MALEG), genes)
pbx <- t(vapply(pts, function(p) Matrix::rowSums(cnt[sg, donor==p, drop=FALSE]), numeric(length(sg)))); colnames(pbx)<-sg
ls <- vapply(pts, function(p) sum(cnt[, donor==p]), numeric(1)); pbx <- log2(sweep(pbx,1,ls,"/")*1e6+1)
zX <- if (FEMALE %in% colnames(pbx)) zscore(pbx[,FEMALE]) else rep(0,length(pts))
yy <- intersect(MALEG,colnames(pbx)); zY <- if(length(yy)) zscore(rowMeans(pbx[,yy,drop=FALSE])) else rep(0,length(pts))
sex <- setNames(ifelse(zY-zX>0.5,"M",ifelse(zY-zX< -0.5,"F","Amb")), pts)

## ---- MAIT pseudobulk per donor (counts) + keep donors with enough MAIT ----
nMAIT <- tapply(is_MAIT, donor, sum)[pts]
keep_d <- pts[!is.na(nMAIT) & nMAIT >= MIN_MAIT & sex[pts] %in% c("M","F")]
cat("Donors with >=", MIN_MAIT, "MAIT and a sex call:", length(keep_d),
    " (M=", sum(sex[keep_d]=="M"), " F=", sum(sex[keep_d]=="F"), ")\n", sep="")
stopifnot(sum(sex[keep_d]=="M")>=3, sum(sex[keep_d]=="F")>=3)
pb <- vapply(keep_d, function(p) Matrix::rowSums(cnt[, donor==p & is_MAIT, drop=FALSE]), numeric(nrow(cnt)))
rownames(pb) <- genes; grp <- factor(sex[keep_d], levels=c("F","M"))

## ---- edgeR quasi-likelihood DE (male vs female) ----
y <- DGEList(pb, group=grp); keepg <- filterByExpr(y, group=grp); y <- y[keepg,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y); des <- model.matrix(~grp); y <- estimateDisp(y, des)
fit <- glmQLFit(y, des); qlf <- glmQLFTest(fit, coef=2)   # coef 2 = M vs F
de <- as.data.frame(topTags(qlf, n=Inf)); de$gene <- rownames(de)
write.csv(de, file.path(OUTDIR,"MAIT_male_vs_female_DE.csv"), row.names=FALSE)
cat("\nTop male-up MAIT genes:", paste(head(de$gene[de$logFC>0][order(de$PValue[de$logFC>0])],10),collapse=", "),"\n")
cat("Top female-up MAIT genes:", paste(head(de$gene[de$logFC<0][order(de$PValue[de$logFC<0])],10),collapse=", "),"\n")
cat("Genes at FDR<0.05:", sum(de$FDR<0.05,na.rm=TRUE),"\n")

## ---- fgsea on ranked genes (Reactome + Hallmark) ----
ranks <- sign(de$logFC) * -log10(pmax(de$PValue,1e-300)); names(ranks) <- de$gene
ranks <- sort(ranks[is.finite(ranks)], decreasing=TRUE)
gs <- function(cat, sub=NULL) { m <- if(is.null(sub)) msigdbr(species="Homo sapiens", category=cat) else
  msigdbr(species="Homo sapiens", category=cat, subcategory=sub)
split(m$gene_symbol, m$gs_name) }
run_fgsea <- function(paths, tag) {
  set.seed(1)                                    # fgsea multilevel is randomised
  fg <- fgsea(paths, ranks, minSize=10, maxSize=500)
  fg <- fg[order(fg$padj), ]; fg$leadingEdge <- vapply(fg$leadingEdge, function(x) paste(head(x,15),collapse=","), character(1))
  write.csv(as.data.frame(fg), file.path(OUTDIR, paste0("MAIT_fgsea_",tag,".csv")), row.names=FALSE)
  cat("\n== fgsea", tag, "– top pathways (padj<0.1) ==\n")
  sig <- fg[fg$padj<0.1, ]; if(nrow(sig)) print(head(sig[,c("pathway","NES","padj")],12)) else cat("(none at padj<0.1)\n")
  fg
}
fg_react <- run_fgsea(gs("C2","CP:REACTOME"), "reactome")
fg_hall  <- run_fgsea(gs("H"),               "hallmark")

## ---- functional-state scores per donor (pseudobulk logCPM) compared by sex ----
logcpm <- edgeR::cpm(y, log=TRUE)   # genes x donors (MAIT pseudobulk)
STATES <- list(
  Cytotoxicity = c("GZMB","GZMA","GZMH","PRF1","GNLY","NKG7","KLRG1","FGFBP2"),
  Cytokine     = c("IFNG","TNF","IL17A","IL17F","IL26","CCL4","CCL5","XCL1"),
  Exhaustion   = c("PDCD1","HAVCR2","LAG3","TIGIT","TOX","ENTPD1","CTLA4"),
  Activation   = c("CD69","HLA-DRA","CD38","MKI67","IL2RA","TNFRSF9"))
sc_state <- sapply(STATES, function(gg){ gg<-gg[gg %in% rownames(logcpm)]
if(!length(gg)) return(rep(NA,ncol(logcpm))); colMeans(t(scale(t(logcpm[gg,,drop=FALSE])))) })
sc_state <- as.data.frame(sc_state); sc_state$sex <- as.character(grp)
st_res <- do.call(rbind, lapply(names(STATES), function(s){
  m<-sc_state[[s]][sc_state$sex=="M"]; f<-sc_state[[s]][sc_state$sex=="F"]
  data.frame(state=s, median_M=round(median(m,na.rm=TRUE),2), median_F=round(median(f,na.rm=TRUE),2),
             higher_in=ifelse(median(m,na.rm=TRUE)>median(f,na.rm=TRUE),"M","F"),
             p=signif(tryCatch(wilcox.test(m,f)$p.value,error=function(e)NA),3), row.names=NULL) }))
write.csv(st_res, file.path(OUTDIR,"MAIT_state_by_sex.csv"), row.names=FALSE)
cat("\n== MAIT functional-state scores by sex ==\n"); print(st_res, row.names=FALSE)

## ---- figures ----
save_fig(function(){                                   # volcano
  par(mar=c(4,4,2,1), cex=0.85)
  plot(de$logFC, -log10(de$PValue), pch=19, cex=.5, col=adjustcolor("grey60",.6),
       xlab="log2 FC (male vs female MAIT)", ylab="-log10 p", main="Male vs female MAIT – DE")
  sigg <- de$FDR<0.05 & !is.na(de$FDR)
  points(de$logFC[sigg], -log10(de$PValue[sigg]), pch=19, cex=.6,
         col=ifelse(de$logFC[sigg]>0,MALE,FEM))
  abline(v=0, col="grey70")
  if(any(sigg)) text(de$logFC[sigg], -log10(de$PValue[sigg]), de$gene[sigg], cex=.5, pos=3)
}, "MAIT_volcano", 6, 5)

topp <- rbind(head(fg_react[fg_react$NES>0,],6), head(fg_react[fg_react$NES<0,],6))
topp <- topp[is.finite(topp$NES),]
if (nrow(topp)) save_fig(function(){
  par(mar=c(4,18,2,1), cex=0.8)
  o <- order(topp$NES); barplot(topp$NES[o], horiz=TRUE, col=ifelse(topp$NES[o]>0,MALE,FEM),
                                names.arg=substr(topp$pathway[o],1,52), las=1, cex.names=0.55, xlab="NES (>0 male-MAIT enriched)",
                                main="Reactome pathways enriched by sex in MAIT")
}, "MAIT_fgsea_reactome", 9, 4.5)

long <- do.call(rbind, lapply(names(STATES), function(s) data.frame(state=s, sex=sc_state$sex, value=sc_state[[s]])))
save_fig(ggplot(long, aes(sex,value,fill=sex)) + geom_boxplot(outlier.shape=NA,alpha=.5) +
           geom_jitter(width=.15,size=1.5,alpha=.8) + facet_wrap(~state,nrow=1,scales="free_y") +
           scale_fill_manual(values=c(M=MALE,F=FEM)) + labs(x=NULL,y="per-donor state score (z)",
                                                            title="MAIT functional-state scores by sex (Andrews healthy)") +
           theme_classic(base_size=9)+theme(legend.position="none"), "SuppFig_S3_MAIT_functional", 9, 3)

cat("\nDone. Wrote DE, fgsea (reactome+hallmark), state tables and figures to", OUTDIR, "\n")