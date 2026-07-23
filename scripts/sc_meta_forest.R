## sc_meta_forest.R -----------------------------------------------------------
## Single-cell meta-analysis of the sex effect across cohorts.
## Reads each cohort's per-donor fraction CSV, computes a common rank-based effect
## size (Cliff's delta, M vs F; >0 = higher in MALES) with bootstrap 95% CI, and
## draws a forest plot with an inverse-variance pooled estimate, for MAIT / Treg /
## CD8. Run on the HPC after the per-cohort validations have written their CSVs.
## ===========================================================================
OUTBASE <- "/data/Blizard-AlazawiLab/rk/MASLD_sex_meta/single_cell"
FIGOUT  <- file.path(OUTBASE, "meta"); dir.create(FIGOUT, showWarnings=FALSE, recursive=TRUE)
MIN_T   <- 50
suppressMessages(library(ggplot2))
set.seed(1)

## cohort registry: name, csv path, optional subset (column=value)
CO <- list(
 list(name="Andrews (healthy)",       csv="andrews2024/out/Andrews2024_per_donor_fractions.csv",     sub=c(group="healthy")),
 list(name="HLiCA (healthy)",         csv="hlica/out/HLiCA_per_patient.csv",                          sub=NULL),
 list(name="Guilliams (lean)",        csv="guilliams/out/Guilliams_clean_per_patient.csv",           sub=c(diet="Lean")),
 list(name="Ramachandran (cirrhotic)",csv="ramachandran/out/Ramachandran2019_per_donor_fractions.csv",sub=c(group="cirrhotic")),
 list(name="Ramachandran (healthy)",  csv="ramachandran/out/Ramachandran2019_per_donor_fractions.csv",sub=c(group="healthy"))
 ## To add your CITE-seq: uncomment and point at a CSV with columns sex, MAIT_pctT, Treg_pctT, CD8_pctT
 ## , list(name="CITE-seq (MASLD)", csv="own_citeseq/out/own_per_donor_fractions.csv", sub=NULL)
)
READOUTS <- c(MAIT_pctT="MAIT", Treg_pctT="Treg", CD8_pctT="CD8")

## Cliff's delta (M vs F) + bootstrap CI + SE
cliffs <- function(m, f) {
 m <- m[!is.na(m)]; f <- f[!is.na(f)]
 if (length(m) < 2 || length(f) < 2) return(c(d=NA, lo=NA, hi=NA, se=NA))
 del <- function(a,b) (sum(outer(a,b,">")) - sum(outer(a,b,"<"))) / (length(a)*length(b))
 d0 <- del(m, f)
 bs <- replicate(2000, del(sample(m, replace=TRUE), sample(f, replace=TRUE)))
 c(d=d0, lo=unname(quantile(bs,.025,na.rm=TRUE)), hi=unname(quantile(bs,.975,na.rm=TRUE)), se=sd(bs,na.rm=TRUE))
}

rows <- list()
for (co in CO) {
 fp <- file.path(OUTBASE, co$csv)
 if (!file.exists(fp)) { cat("missing:", co$csv, "\n"); next }
 x <- read.csv(fp, stringsAsFactors=FALSE)
 names(x)[names(x)=="n_Tcells"] <- "n_T"                 # standardise capture column
 if (!is.null(co$sub)) x <- x[as.character(x[[names(co$sub)]]) == unname(co$sub), , drop=FALSE]
 if ("n_T" %in% names(x)) x <- x[!is.na(x$n_T) & x$n_T >= MIN_T, ]
 x <- x[x$sex %in% c("M","F"), ]
 nM <- sum(x$sex=="M"); nF <- sum(x$sex=="F")
 for (rv in names(READOUTS)) {
  if (!rv %in% names(x)) next
  ef <- cliffs(x[[rv]][x$sex=="M"], x[[rv]][x$sex=="F"])
  p  <- tryCatch(wilcox.test(x[[rv]][x$sex=="M"], x[[rv]][x$sex=="F"])$p.value, error=function(e) NA)
  rows[[length(rows)+1]] <- data.frame(readout=READOUTS[[rv]], cohort=co$name,
                                       nM=nM, nF=nF, delta=ef["d"], lo=ef["lo"], hi=ef["hi"], se=ef["se"], p=p, row.names=NULL)
 }
}
res <- do.call(rbind, rows)

## inverse-variance pooled delta per readout (random-effects light)
pooled <- do.call(rbind, lapply(unique(res$readout), function(rd) {
 s <- res[res$readout==rd & is.finite(res$delta) & is.finite(res$se) & res$se>0, ]
 if (!nrow(s)) return(NULL)
 w <- 1/s$se^2; d <- sum(w*s$delta)/sum(w); se <- sqrt(1/sum(w))
 data.frame(readout=rd, cohort="POOLED", nM=sum(s$nM), nF=sum(s$nF),
            delta=d, lo=d-1.96*se, hi=d+1.96*se, se=se, p=2*pnorm(-abs(d/se)), row.names=NULL)
}))
res <- rbind(res, pooled)
write.csv(res, file.path(FIGOUT, "sc_meta_effects.csv"), row.names=FALSE)
cat("\n==== single-cell meta-analysis (Cliff's delta, >0 = higher in MALES) ====\n")
print(res[order(res$readout, res$cohort=="POOLED"), c("readout","cohort","nM","nF","delta","lo","hi","p")], row.names=FALSE, digits=2)

## forest plot
res$readout <- factor(res$readout, levels=c("MAIT","Treg","CD8"))
res$cohort  <- factor(res$cohort, levels=rev(c(setdiff(unique(res$cohort),"POOLED"),"POOLED")))
res$col <- ifelse(res$cohort=="POOLED","#111111", ifelse(res$delta>0,"#2c7fb8","#d95f0e"))
p <- ggplot(res, aes(delta, cohort)) +
 geom_vline(xintercept=0, linewidth=.4, colour="grey60") +
 geom_errorbarh(aes(xmin=lo, xmax=hi), height=.22, colour=res$col) +
 geom_point(aes(size=ifelse(cohort=="POOLED",3.2,2.2)), colour=res$col, shape=ifelse(res$cohort=="POOLED",18,16)) +
 facet_wrap(~readout, nrow=1) + scale_size_identity() +
 coord_cartesian(xlim=c(-1,1)) +
 labs(x="Cliff's delta  (<0 higher in females ← | → higher in males >0)", y=NULL,
      title="Single-cell meta-analysis of the hepatic immune sex effect") +
 theme_classic(base_size=10) + theme(panel.grid.major.y=element_line(colour="grey92"))
ggsave(file.path(FIGOUT, "sc_meta_forest.png"), p, width=10, height=3.4, dpi=300)
cat("\nWrote", file.path(FIGOUT,"sc_meta_forest.png"), "and sc_meta_effects.csv\n")