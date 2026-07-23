## sc_guilliams_clean.R -------------------------------------------------------
## Corrected sex analysis of the Guilliams human Liver Cell Atlas.
## Fixes the first-pass confounds:
##   - donor = PATIENT (not sample)  -> no pseudoreplication
##   - single-CELL only (drop single-nucleus)
##   - cell types from the atlas 'annot' (cDC = cDC1s/cDC2s/Mig.cDCs); MAIT/Treg/
##     CD8 gated by markers WITHIN annot=="T cells"
##   - require >=50 T cells / patient so no fraction rests on a few cells
##   - stratify by diet (Lean vs Obese) and expose the sex x diet confound
## ===========================================================================
OBJ_PATH <- "/data/Blizard-AlazawiLab/rk/MASLD_sex_meta/single_cell/guilliams/GSE192742_seurat.rds"
ANNOT    <- "/gpfs/scratch/hdx044/guilliams/annot_humanAll.csv"
OUTDIR   <- "/data/Blizard-AlazawiLab/rk/MASLD_sex_meta/single_cell/guilliams/out"
MIN_T    <- 50    # minimum T cells per patient to trust a fraction
## ===========================================================================
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
suppressMessages({library(Seurat); library(Matrix); library(data.table)})

obj <- readRDS(OBJ_PATH); DefaultAssay(obj) <- "RNA"
an  <- as.data.frame(fread(ANNOT)); rownames(an) <- as.character(an$cell)
cells <- colnames(obj); an <- an[cells, ]
obj$patient    <- an$patient
obj$diet       <- an$diet
obj$annot      <- as.character(an$annot)
obj$typeSample <- an$typeSample

## keep single-CELL only (drop single-nucleus, which can't see MAIT)
obj <- obj[, obj$typeSample %in% c("scRnaSeq","citeSeq") & !is.na(obj$patient)]
cnt <- tryCatch(GetAssayData(obj, assay="RNA", layer="counts"), error=function(e) GetAssayData(obj, assay="RNA", slot="counts"))
patient <- obj$patient; diet <- obj$diet; annot <- obj$annot
pts <- sort(unique(patient))
cat("Single-cell patients:", length(pts), "| cells:", ncol(obj), "\n")

genes <- rownames(cnt); has <- function(g) g %in% genes
gvec <- function(g) if (has(g)) as.numeric(cnt[g, ]) else rep(0, ncol(cnt))
anyMark <- function(gs){ gs<-gs[vapply(gs,has,logical(1))]; if(!length(gs)) return(rep(FALSE,ncol(cnt))); Reduce(`|`, lapply(gs, function(g) gvec(g)>0)) }

## cell-type masks: annot for lineages, markers within T cells for subsets
is_imm  <- !annot %in% c("Hepatocytes","Endothelial cells","Fibroblasts","Cholangiocytes")
is_T    <- annot == "T cells"
is_MAIT <- is_T & anyMark(c("SLC4A10","TRAV1-2"))
is_CD8  <- is_T & anyMark(c("CD8A","CD8B"))
is_Treg <- is_T & anyMark(c("FOXP3"))
is_cDC  <- annot %in% c("cDC1s","cDC2s","Mig.cDCs")

## sex per patient (pseudobulk XIST vs Y)
zscore <- function(x){ m<-mean(x,na.rm=TRUE); s<-sqrt(mean((x-m)^2,na.rm=TRUE)); if(is.na(s)||s==0) x*0 else (x-m)/s }
FEMALE<-"XIST"; MALE<-c("RPS4Y1","DDX3Y","EIF1AY","UTY","KDM5D","USP9Y","NLGN4Y","ZFY","TXLNGY")
sg <- intersect(c(FEMALE,MALE), genes)
pb <- t(vapply(pts, function(p) Matrix::rowSums(cnt[sg, patient==p, drop=FALSE]), numeric(length(sg)))); colnames(pb)<-sg
ls <- vapply(pts, function(p) sum(cnt[, patient==p]), numeric(1)); pb <- log2(sweep(pb,1,ls,"/")*1e6+1)
zX <- if (FEMALE %in% colnames(pb)) zscore(pb[,FEMALE]) else rep(0,length(pts))
yy <- intersect(MALE,colnames(pb)); zY <- if(length(yy)) zscore(rowMeans(pb[,yy,drop=FALSE])) else rep(0,length(pts))
sex <- ifelse(zY-zX>0.5,"M",ifelse(zY-zX< -0.5,"F","Amb")); names(sex)<-pts

frac <- function(num,den){ n<-tapply(num,patient,sum)[pts]; d<-tapply(den,patient,sum)[pts]; ifelse(d>0,100*n/d,NA) }
df <- data.frame(patient=pts, sex=sex[pts],
                 diet=tapply(diet,patient,function(x) x[1])[pts],
                 n_T=as.integer(tapply(is_T,patient,sum)[pts]), n_MAIT=as.integer(tapply(is_MAIT,patient,sum)[pts]),
                 MAIT_pctT=frac(is_MAIT,is_T), Treg_pctT=frac(is_Treg,is_T), CD8_pctT=frac(is_CD8,is_T), cDC_pctImm=frac(is_cDC,is_imm),
                 row.names=NULL)
df <- df[df$sex %in% c("M","F"), ]

cat("\n--- sex x diet (patient level) ---\n"); print(table(df$sex, df$diet))
cat("\n--- T-cell capture per patient (median) by sex ---\n"); print(tapply(df$n_T, df$sex, median))
usable <- df[df$n_T >= MIN_T, ]
cat("\nPatients with >=", MIN_T, "T cells:", nrow(usable), "/", nrow(df),
    " (M=", sum(usable$sex=="M"), " F=", sum(usable$sex=="F"), ")\n", sep="")
write.csv(df, file.path(OUTDIR,"Guilliams_clean_per_patient.csv"), row.names=FALSE)

EXPECT <- c(MAIT_pctT="M", Treg_pctT="F", CD8_pctT="F", cDC_pctImm="F")
report <- function(d2,label){ d2<-d2[d2$sex %in% c("M","F"),]
if (sum(d2$sex=="M")<2 || sum(d2$sex=="F")<2){cat("\n[",label,"] too few per sex (M=",sum(d2$sex=="M")," F=",sum(d2$sex=="F"),") - skipped\n",sep="");return(invisible())}
res<-do.call(rbind, lapply(names(EXPECT), function(v){ m<-d2[[v]][d2$sex=="M"]; f<-d2[[v]][d2$sex=="F"]
p<-tryCatch(wilcox.test(m,f)$p.value,error=function(e) NA); hi<-if(median(m,na.rm=TRUE)>median(f,na.rm=TRUE))"M" else "F"
data.frame(readout=v, med_M=round(median(m,na.rm=TRUE),2), med_F=round(median(f,na.rm=TRUE),2),
           higher_in=hi, predicted=EXPECT[[v]], matches=ifelse(hi==EXPECT[[v]],"YES","no"), p=signif(p,3), row.names=NULL)}))
cat("\n==== ",label," (M=",sum(d2$sex=="M"),", F=",sum(d2$sex=="F")," patients) ====\n",sep="")
print(res,row.names=FALSE) }
report(usable, "all usable single-cell patients")
for (dg in unique(na.omit(usable$diet))) report(usable[usable$diet==dg, ], paste("diet:", dg))
cat("\nDone. If capture is still badly sex-imbalanced or a diet stratum is empty, Guilliams\n")
cat("cannot give a clean sex contrast and should be set aside with that reason stated.\n")