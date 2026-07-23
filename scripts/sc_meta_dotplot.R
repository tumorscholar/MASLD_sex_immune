## sc_meta_dotplot.R ----------------------------------------------------------
## Combined single-cell validation figure (manuscript Figure 7): per-donor DOT
## PLOTS by sex, one panel per readout (rows) x cohort (columns). The first
## column is our in-house CITE-seq cohort (section 3.7); the remaining columns
## are the four public cohorts. Each dot = one donor, black crossbar = median,
## blue = male / orange = female. Also exports the pooled per-donor table.
##
## Reads the per-cohort fraction CSVs written by the validations, plus the
## in-house fractions (data/own_cohort_percell_fractions.csv). Column names are
## harmonised: proportion columns (_ofT/_ofImm) or percent columns
## (_pctT/_pctImm) are both accepted and shown as percentages.
## ===========================================================================
OUTBASE <- Sys.getenv("MASLD_SC",
             "/path/to/MASLD_sex_meta/single_cell")
FIGOUT  <- file.path(OUTBASE, "meta"); dir.create(FIGOUT, showWarnings=FALSE, recursive=TRUE)
## where the in-house fractions live (shipped in the repo under data/)
OWN_CSV <- Sys.getenv("OWN_CITESEQ_CSV", "")
if (OWN_CSV == "") {
  for (p in c(file.path(OUTBASE, "own_citeseq/own_cohort_percell_fractions.csv"),
              "../data/own_cohort_percell_fractions.csv",
              "data/own_cohort_percell_fractions.csv"))
    if (file.exists(p)) { OWN_CSV <- p; break }
}
MIN_T <- 50
suppressMessages(library(ggplot2))
MALE_COL <- "#2c7fb8"; FEM_COL <- "#d95f0e"

## cohort list: in-house first, then the four public cohorts.
## `path` is absolute for the in-house file, else relative to OUTBASE.
CO <- list(
  list(name="In-house\nCITE-seq",      path=OWN_CSV,                                                     sub=NULL, abs=TRUE),
  list(name="Andrews\n(healthy)",       path="andrews2024/out/Andrews2024_per_donor_fractions.csv",      sub=c(group="healthy")),
  list(name="HLiCA\n(healthy)",         path="hlica/out/HLiCA_per_patient.csv",                           sub=NULL),
  list(name="Guilliams\n(lean)",        path="guilliams/out/Guilliams_clean_per_patient.csv",            sub=c(diet="Lean")),
  list(name="Ramachandran\n(cirrhotic)",path="ramachandran/out/Ramachandran2019_per_donor_fractions.csv",sub=c(group="cirrhotic"))
)
## readout label -> accepted column names (percent first, then proportion)
READOUTS <- list(
  "MAIT (% T)"   = c("MAIT_pctT","MAIT_ofT"),
  "Treg (% T)"   = c("Treg_pctT","Treg_ofT"),
  "CD8 (% T)"    = c("CD8_pctT","CD8_ofT"),
  "cDC (% imm)"  = c("cDC_pctImm","cDC_ofImm"))

long <- list(); stats <- list()
for (co in CO) {
  fp <- if (isTRUE(co$abs)) co$path else file.path(OUTBASE, co$path)
  if (is.null(fp) || fp == "" || !file.exists(fp)) { cat("missing:", co$name, "->", fp, "\n"); next }
  x <- read.csv(fp, stringsAsFactors=FALSE)
  names(x)[names(x)=="n_Tcells"] <- "n_T"
  if (!is.null(co$sub) && names(co$sub) %in% names(x))
    x <- x[as.character(x[[names(co$sub)]]) == unname(co$sub), , drop=FALSE]
  if ("n_T" %in% names(x)) x <- x[!is.na(x$n_T) & x$n_T >= MIN_T, ]
  x <- x[x$sex %in% c("M","F"), ]
  if (sum(x$sex=="M") < 2 || sum(x$sex=="F") < 2) { cat("skip (too few/sex):", co$name, "\n"); next }
  for (lab in names(READOUTS)) {
    col <- READOUTS[[lab]][READOUTS[[lab]] %in% names(x)][1]
    if (is.na(col)) next
    v <- x[[col]]; if (max(v, na.rm=TRUE) <= 1.5) v <- v * 100   # proportion -> %
    long[[length(long)+1]]  <- data.frame(cohort=co$name, readout=lab, sex=x$sex, value=v)
    pv <- tryCatch(wilcox.test(v[x$sex=="M"], v[x$sex=="F"])$p.value, error=function(e) NA)
    stats[[length(stats)+1]] <- data.frame(cohort=co$name, readout=lab, p=pv)
  }
}
D <- do.call(rbind, long); S <- do.call(rbind, stats)
D <- D[!is.na(D$value), ]
write.csv(D, file.path(FIGOUT, "sc_meta_donor_data.csv"), row.names=FALSE)

lev_co <- vapply(CO, function(z) z$name, character(1))
lev_co <- lev_co[lev_co %in% D$cohort]
D$cohort  <- factor(D$cohort,  levels=lev_co);  S$cohort  <- factor(S$cohort,  levels=lev_co)
D$readout <- factor(D$readout, levels=names(READOUTS)); S$readout <- factor(S$readout, levels=names(READOUTS))
S$lab <- ifelse(is.na(S$p), "", paste0("p=", signif(S$p,2)))
ytop <- aggregate(value~cohort+readout, D, function(v) max(v)*1.05)
S <- merge(S, ytop, all.x=TRUE)

p <- ggplot(D, aes(sex, value, colour=sex)) +
  geom_jitter(width=0.18, height=0, size=1.5, alpha=0.75) +
  stat_summary(fun=median, geom="crossbar", width=0.55, linewidth=0.3, colour="black") +
  geom_text(data=S, aes(x=1.5, y=value, label=lab), inherit.aes=FALSE, size=2.4, vjust=0, colour="grey25") +
  facet_grid(readout ~ cohort, scales="free_y", switch="y") +
  scale_colour_manual(values=c(M=MALE_COL, F=FEM_COL)) +
  labs(x=NULL, y="per-donor fraction (%)",
       title="Single-cell validation of the hepatic immune sex effect: in-house and public cohorts",
       subtitle="each dot = one donor · black bar = median · blue male / orange female") +
  theme_bw(base_size=9) +
  theme(legend.position="none", panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="grey95", colour=NA),
        strip.placement="outside", plot.subtitle=element_text(size=8, colour="grey40"))

## save the combined Figure 7 alongside Figs 1-6 in the main figures folder,
## in submission formats: PDF (vector) + TIFF (300 dpi, LZW) + PNG
FIGDIR <- Sys.getenv("MASLD_FIGDIR", file.path(dirname(OUTBASE), "figures"))
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
NM <- "Fig7_singlecell_combined"
ggsave(file.path(FIGDIR, paste0(NM,".pdf")),  p, width=11, height=6.5)
ggsave(file.path(FIGDIR, paste0(NM,".tiff")), p, width=11, height=6.5, dpi=300, compression="lzw")
ggsave(file.path(FIGDIR, paste0(NM,".png")),  p, width=11, height=6.5, dpi=300)
cat("Wrote", NM, "(pdf, tiff, png) to", FIGDIR, "and sc_meta_donor_data.csv to", FIGOUT, "\n")
cat("Columns plotted (left to right):", paste(levels(D$cohort), collapse=" | "), "\n")
