## sc_build_andrews.R ---------------------------------------------------------
## Build one Seurat object from the Andrews 2024 (GSE243981) SINGLE-CELL libraries.
## The data are matrix-market triplets (barcodes/genes/matrix.mtx.gz), loaded with
## ReadMtx. We keep only _SC_ (single-cell) libraries; single-nucleus (_SN_) and
## Visium (spatial) samples are skipped. Each donor is tagged healthy / PSC / PBC.
## No clustering/annotation needed – the validation gates cell types from markers.
## Run in RStudio on the HPC after untarring GSE243981_RAW.tar.
## ===========================================================================
IN_DIR   <- "/gpfs/scratch/hdx044/andrews2024"   # RAW data on SCRATCH (extracted tar)
OUT_RDS  <- "/data/Blizard-AlazawiLab/rk/MASLD_sex_meta/single_cell/andrews2024/GSE243981_seurat.rds"
MIN_GENES <- 200; MAX_GENES <- 6000; MAX_MT <- 20
## ===========================================================================
suppressMessages({library(Seurat); library(Matrix)})
dir.create(dirname(OUT_RDS), showWarnings = FALSE, recursive = TRUE)

files <- list.files(IN_DIR, full.names = TRUE, recursive = TRUE)
bn <- basename(files)
## single-cell matrix files only (exclude single-nucleus, Visium, and any .h5)
mtx <- files[grepl("_SC_", bn) & grepl("matrix\\.mtx\\.gz$", bn) &
              !grepl("_SN_|VISIUM|\\.h5$", bn)]
cat("Single-cell libraries found:", length(mtx), "\n")
stopifnot(length(mtx) >= 1)

objs <- list()
for (mp in mtx) {
 gsm <- sub("_.*", "", basename(mp))                       # e.g. GSM7802679
 grp <- files[startsWith(basename(files), paste0(gsm, "_"))]
 cells <- grp[grepl("barcodes\\.tsv\\.gz$", basename(grp))]
 feat  <- grp[grepl("(genes|features)\\.tsv\\.gz$", basename(grp))]
 if (!length(cells) || !length(feat)) { cat("  skip", gsm, "(missing barcodes/features)\n"); next }
 sid   <- sub("_barcodes.*", "", sub(paste0("^", gsm, "_"), "", basename(cells[1])))  # C41_NPC_SC_3pr
 donor <- sub("_(NPC_)?(flush_)?SC_.*", "", sid)           # C41 / C51 / PSC012 ...
 group <- if (grepl("^C", donor)) "healthy" else if (grepl("^PSC", donor)) "PSC" else if (grepl("^PBC", donor)) "PBC" else "other"
 m <- tryCatch(ReadMtx(mtx = mp, cells = cells[1], features = feat[1], feature.column = 2),
               error = function(e) tryCatch(ReadMtx(mtx = mp, cells = cells[1], features = feat[1], feature.column = 1),
                                            error = function(e2) { cat("  read fail", sid, ":", conditionMessage(e2), "\n"); NULL }))
 if (is.null(m)) next
 o <- CreateSeuratObject(m, project = sid, min.cells = 3, min.features = MIN_GENES)
 o$donor <- donor; o$group <- group; o$library <- sid
 o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
 o <- subset(o, subset = nFeature_RNA >= MIN_GENES & nFeature_RNA <= MAX_GENES & percent.mt <= MAX_MT)
 objs[[sid]] <- o
 cat(sprintf("  %-22s donor=%-8s group=%-8s %6d cells\n", sid, donor, group, ncol(o)))
}
stopifnot(length(objs) >= 1)
obj <- if (length(objs) == 1) objs[[1]] else
 merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), project = "Andrews2024")
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)   # Seurat v5: one counts layer

cat("\nMerged:", ncol(obj), "cells |", length(unique(obj$donor)), "donors\n")
cat("Donors by group:\n"); print(tapply(obj$donor, obj$group, function(x) length(unique(x))))
saveRDS(obj, OUT_RDS)
cat("\nWrote", OUT_RDS, "\n")
cat("Next: in sc_validate_andrews.R set DONOR_COL='donor', GROUP_COL='group', SEX_COL=NA, then run.\n")
