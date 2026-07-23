## sc_build_ramachandran.R ----------------------------------------------------
## Build one Seurat object from Ramachandran et al. 2019 (GSE136103), keeping only
## the HUMAN LIVER, CD45+ (immune-enriched) libraries. Drops CD45- fractions,
## blood/PBMC, and mouse samples. Each donor tagged healthy / cirrhotic.
## Data are matrix-market triplets (barcodes/genes/matrix.mtx.gz) -> ReadMtx.
## Download first:
##   mkdir -p /path/to/scratch/ramachandran
##   cd /path/to/scratch/ramachandran
##   wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/suppl/GSE136103_RAW.tar'
##   tar -xvf GSE136103_RAW.tar
## ===========================================================================
IN_DIR   <- file.path(Sys.getenv("MASLD_SCRATCH", "/path/to/scratch"), "ramachandran")
OUT_RDS  <- file.path(Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta"), "single_cell/ramachandran/GSE136103_seurat.rds")
MIN_GENES <- 200; MAX_GENES <- 6000; MAX_MT <- 20
## ===========================================================================
suppressMessages({library(Seurat); library(Matrix)})
dir.create(dirname(OUT_RDS), showWarnings = FALSE, recursive = TRUE)

files <- list.files(IN_DIR, full.names = TRUE, recursive = TRUE)
bn <- basename(files)
mtx <- files[grepl("matrix\\.mtx\\.gz$", bn)]
cat("Matrix files in archive:", length(mtx), "\n")

## keep human liver CD45+ only: title has healthy/cirrhotic + cd45 (positive),
## and is NOT a cd45- fraction, blood, or mouse sample
keep_sample <- function(x) grepl("healthy|cirrhotic", x, ignore.case=TRUE) &
 grepl("cd45", x, ignore.case=TRUE) &
 !grepl("cd45-|cd45neg|blood|mouse|pbmc", x, ignore.case=TRUE)

objs <- list()
for (mp in mtx) {
 gsm <- sub("_.*", "", basename(mp))
 grp <- files[startsWith(basename(files), paste0(gsm, "_"))]
 cells <- grp[grepl("barcodes\\.tsv\\.gz$", basename(grp))]
 feat  <- grp[grepl("(genes|features)\\.tsv\\.gz$", basename(grp))]
 if (!length(cells) || !length(feat)) next
 sid <- sub("_barcodes.*", "", sub(paste0("^", gsm, "_"), "", basename(cells[1])))  # e.g. Healthy1_Cd45+
 if (!keep_sample(sid)) next
 donor <- sub("_[Cc][Dd]45.*", "", sid)                                              # Healthy1 / Cirrhotic3
 group <- if (grepl("^healthy", donor, ignore.case=TRUE)) "healthy" else "cirrhotic"
 m <- tryCatch(ReadMtx(mtx=mp, cells=cells[1], features=feat[1], feature.column=2),
               error=function(e) tryCatch(ReadMtx(mtx=mp, cells=cells[1], features=feat[1], feature.column=1),
                                          error=function(e2){cat("  read fail",sid,":",conditionMessage(e2),"\n");NULL}))
 if (is.null(m)) next
 o <- CreateSeuratObject(m, project=sid, min.cells=3, min.features=MIN_GENES)
 o$donor <- donor; o$group <- group; o$library <- sid
 o[["percent.mt"]] <- PercentageFeatureSet(o, pattern="^MT-")
 o <- subset(o, subset = nFeature_RNA>=MIN_GENES & nFeature_RNA<=MAX_GENES & percent.mt<=MAX_MT)
 objs[[sid]] <- o
 cat(sprintf("  %-22s donor=%-12s group=%-10s %6d cells\n", sid, donor, group, ncol(o)))
}
stopifnot(length(objs) >= 1)
obj <- if (length(objs)==1) objs[[1]] else merge(objs[[1]], y=objs[-1], add.cell.ids=names(objs), project="Ramachandran2019")
obj <- tryCatch(JoinLayers(obj), error=function(e) obj)
cat("\nMerged:", ncol(obj), "cells |", length(unique(obj$donor)), "donors\n")
cat("Donors by group:\n"); print(tapply(obj$donor, obj$group, function(x) length(unique(x))))
saveRDS(obj, OUT_RDS, compress = FALSE)   # compress=FALSE = much faster save
cat("\nWrote", OUT_RDS, "\n")
cat("Next: run sc_validate_ramachandran.R\n")
