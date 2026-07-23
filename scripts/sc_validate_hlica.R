## sc_validate_hlica.R --------------------------------------------------------
## Sex validation on HLiCA – the integrated HEALTHY human liver atlas (CELLxGENE).
## Healthy + both sexes => a balanced reinforcement of the CONSTITUTIONAL male-MAIT
## prediction. Uses CELLxGENE's standardized metadata (donor_id, sex, cell_type),
## so sex is RECORDED (we also cross-check it from expression).
##
## DOWNLOAD (to scratch): on cellxgene.cziscience.com, open the HLiCA / integrated
## healthy human liver atlas dataset -> Download -> choose the Seurat (.rds) file
## (or .h5ad). Put it in /gpfs/scratch/hdx044/hlica/ and point OBJ_PATH at it.
## ===========================================================================
OBJ_PATH <- "/gpfs/scratch/hdx044/hlica/hlica_lymphocyte.h5ad"   # CELLxGENE .h5ad (or .rds)
OUTDIR   <- "/data/Blizard-AlazawiLab/rk/MASLD_sex_meta/single_cell/hlica/out"
MIN_T    <- 50
## ===========================================================================
DATASET <- "HLiCA"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
suppressMessages({library(Seurat); library(Matrix)})

## --- load (.rds direct; .h5ad via schard, pure-R) ---
if (grepl("\\.h5ad$", OBJ_PATH)) {
 if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
 if (!requireNamespace("schard",  quietly=TRUE)) remotes::install_github("cellgeni/schard", upgrade="never")
 obj <- tryCatch(schard::h5ad2seurat(OBJ_PATH, use.raw = TRUE),   # use.raw = raw counts
                 error = function(e) schard::h5ad2seurat(OBJ_PATH))
} else obj <- readRDS(OBJ_PATH)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data
cat("Cells:", ncol(obj), "| metadata columns:", paste(head(names(md),40), collapse=", "), "\n")

## CELLxGENE standard fields (with fallbacks)
pick <- function(cands) { h <- names(md)[tolower(names(md)) %in% cands]; if (length(h)) h[1] else NA }
donorc <- pick(c("donor_id","donor","patient","sampleid","sample"))
sexc   <- pick(c("sex","gender"))
ctc    <- pick(c("cell_type","celltype","annotation","annot","cell_type_annotation"))
stopifnot("no donor column found" = !is.na(donorc))
cat("donor:", donorc, "| sex:", sexc, "| cell_type:", ctc, "\n")
donor <- as.character(md[[donorc]]); pts <- sort(unique(donor))
ct    <- if (!is.na(ctc)) tolower(as.character(md[[ctc]])) else rep("", nrow(md))

cnt <- tryCatch(GetAssayData(obj, assay="RNA", layer="counts"), error=function(e) GetAssayData(obj, assay="RNA", slot="counts"))
## CELLxGENE h5ad indexes genes by Ensembl ID -> remap rownames to gene symbols
if (any(grepl("^ENSG", head(rownames(cnt), 20)))) {
 fm <- tryCatch(obj[["RNA"]][[]], error=function(e) tryCatch(GetAssay(obj,"RNA")@meta.features, error=function(e2) NULL))
 symcol <- if (!is.null(fm)) intersect(c("feature_name","gene_symbols","gene_symbol","symbol","gene_name","Gene"), colnames(fm)) else character(0)
 if (length(symcol)) {
  sym <- as.character(fm[rownames(cnt), symcol[1]])
  keep <- !is.na(sym) & sym != "" & !duplicated(sym)
  cnt <- cnt[keep, ]; rownames(cnt) <- sym[keep]
  cat("Remapped", sum(keep), "Ensembl IDs to gene symbols via '", symcol[1], "'\n", sep="")
 } else cat("WARNING: rownames look like Ensembl but no symbol column found – markers may not match.\n")
}
genes <- rownames(cnt); has <- function(g) g %in% genes
gv <- function(g) if (has(g)) as.numeric(cnt[g,]) else rep(0, ncol(cnt))
anyMark <- function(gs){ gs<-gs[vapply(gs,has,logical(1))]; if(!length(gs)) return(rep(FALSE,ncol(cnt))); Reduce(`|`, lapply(gs,function(g) gv(g)>0)) }

## cell-type masks: prefer cell_type labels, fall back to markers
is_T   <- grepl("t cell|t-cell|cd4|cd8|regulatory t|mait", ct) | anyMark(c("CD3D","CD3E","CD3G"))
is_imm <- !grepl("hepatocyte|endothel|stellate|cholangio|fibroblast|mesenchym", ct) & (grepl("cell", ct) | is_T)
is_MAIT<- (grepl("mucosal.*invariant|mait", ct)) | (is_T & anyMark(c("SLC4A10","TRAV1-2")))
is_CD8 <- (grepl("cd8", ct)) | (is_T & anyMark(c("CD8A","CD8B")))
is_Treg<- (grepl("regulatory t", ct)) | (is_T & anyMark(c("FOXP3")))
is_cDC <- grepl("conventional dendritic|cdc|myeloid dendritic", ct) | (grepl("dendritic", ct) & !grepl("plasmacytoid|pdc", ct)) |
 (anyMark(c("CD1C","CLEC9A","CLEC10A")) & !anyMark(c("LILRA4","IL3RA")))

## sex: recorded, cross-checked against expression
zscore <- function(x){ m<-mean(x,na.rm=TRUE); s<-sqrt(mean((x-m)^2,na.rm=TRUE)); if(is.na(s)||s==0) x*0 else (x-m)/s }
FEMALE<-"XIST"; MALE<-c("RPS4Y1","DDX3Y","EIF1AY","UTY","KDM5D","USP9Y","NLGN4Y","ZFY","TXLNGY")
sg<-intersect(c(FEMALE,MALE),genes)
pb<-t(vapply(pts,function(p) Matrix::rowSums(cnt[sg,donor==p,drop=FALSE]),numeric(length(sg)))); colnames(pb)<-sg
ls<-vapply(pts,function(p) sum(cnt[,donor==p]),numeric(1)); pb<-log2(sweep(pb,1,ls,"/")*1e6+1)
zX<-if(FEMALE%in%colnames(pb)) zscore(pb[,FEMALE]) else rep(0,length(pts))
yy<-intersect(MALE,colnames(pb)); zY<-if(length(yy)) zscore(rowMeans(pb[,yy,drop=FALSE])) else rep(0,length(pts))
sex_expr<-ifelse(zY-zX>0.5,"M",ifelse(zY-zX< -0.5,"F","Amb")); names(sex_expr)<-pts
if (!is.na(sexc)) {
 rec<-tapply(tolower(as.character(md[[sexc]])),donor,function(x) x[1])[pts]
 rec<-ifelse(grepl("^m|male",rec) & !grepl("female",rec),"M",ifelse(grepl("^f|female",rec),"F",NA))
 ok<-!is.na(rec)&sex_expr!="Amb"; if(any(ok)) cat(sprintf("Sex concordance (recorded vs expression): %.1f%% (%d)\n",100*mean(rec[ok]==sex_expr[ok]),sum(ok)))
 sex<-ifelse(!is.na(rec),rec,sex_expr)   # trust recorded, fill gaps from expression
} else sex<-sex_expr
names(sex)<-pts

frac<-function(num,den){ n<-tapply(num,donor,sum)[pts]; d<-tapply(den,donor,sum)[pts]; ifelse(d>0,100*n/d,NA) }
df<-data.frame(patient=pts, sex=sex[pts], n_T=as.integer(tapply(is_T,donor,sum)[pts]), n_MAIT=as.integer(tapply(is_MAIT,donor,sum)[pts]),
               MAIT_pctT=frac(is_MAIT,is_T), Treg_pctT=frac(is_Treg,is_T), CD8_pctT=frac(is_CD8,is_T), cDC_pctImm=frac(is_cDC,is_imm), row.names=NULL)
df<-df[df$sex %in% c("M","F"),]
write.csv(df, file.path(OUTDIR,"HLiCA_per_patient.csv"), row.names=FALSE)
cat("\nPatients:", nrow(df), "(M=",sum(df$sex=="M")," F=",sum(df$sex=="F"),") | MAIT median:",median(df$n_MAIT,na.rm=TRUE),"\n")
usable<-df[df$n_T>=MIN_T,]; cat("With >=",MIN_T,"T cells:",nrow(usable)," (M=",sum(usable$sex=="M")," F=",sum(usable$sex=="F"),")\n",sep="")

EXPECT<-c(MAIT_pctT="M",Treg_pctT="F",CD8_pctT="F",cDC_pctImm="F")
d2<-usable[usable$sex%in%c("M","F"),]
if (sum(d2$sex=="M")>=2 && sum(d2$sex=="F")>=2) {
 res<-do.call(rbind,lapply(names(EXPECT),function(v){ m<-d2[[v]][d2$sex=="M"]; f<-d2[[v]][d2$sex=="F"]
 p<-tryCatch(wilcox.test(m,f)$p.value,error=function(e) NA); hi<-if(median(m,na.rm=TRUE)>median(f,na.rm=TRUE))"M" else "F"
 data.frame(readout=v,med_M=round(median(m,na.rm=TRUE),2),med_F=round(median(f,na.rm=TRUE),2),higher_in=hi,predicted=EXPECT[[v]],matches=ifelse(hi==EXPECT[[v]],"YES","no"),p=signif(p,3),row.names=NULL)}))
 cat("\n==== ",DATASET," (M=",sum(d2$sex=="M"),", F=",sum(d2$sex=="F")," patients) ====\n",sep=""); print(res,row.names=FALSE)
 cat("Directions reproduced:",sum(res$matches=="YES"),"/ 4\n")
 write.csv(res, file.path(OUTDIR,"HLiCA_sex_validation.csv"), row.names=FALSE)
} else cat("\nToo few per sex after capture filter – report counts only.\n")
cat("\nWrote to", OUTDIR, "\n")
