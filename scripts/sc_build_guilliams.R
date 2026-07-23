## sc_build_guilliams.R -------------------------------------------------------
## Build a Seurat object from the Guilliams 2022 Liver Cell Atlas HUMAN data.
## Unlike the GEO per-sample tarballs, this is ONE matrix + an annotation CSV.
## Download first (to scratch):
##   mkdir -p /path/to/scratch/guilliams && cd /path/to/scratch/guilliams
##   wget https://www.livercellatlas.org/data_files/toDownload/rawData_human.zip
##   wget https://www.livercellatlas.org/data_files/toDownload/annot_humanAll.csv
##   unzip rawData_human.zip
## Then run this. It auto-detects the matrix files and the donor/group columns;
## check the "Annotation columns" print and adjust DONOR_COL/GROUP_COL if needed.
## ===========================================================================
IN_DIR   <- file.path(Sys.getenv("MASLD_SCRATCH", "/path/to/scratch"), "guilliams")
ANNOT    <- file.path(IN_DIR, "annot_humanAll.csv")
OUT_RDS  <- file.path(Sys.getenv("MASLD_REALDIR", "/path/to/MASLD_sex_meta"), "single_cell/guilliams/GSE192742_seurat.rds")
DONOR_COL <- NA   # set to a column name after you see the print; NA = auto-detect (patient/sample/…)
GROUP_COL <- NA   # disease/diet/condition column if present; NA = auto-detect, else "all"
MIN_GENES <- 200; MAX_GENES <- 6000; MAX_MT <- 20
## ===========================================================================
suppressMessages({library(Seurat); library(Matrix); library(data.table)})
dir.create(dirname(OUT_RDS), showWarnings = FALSE, recursive = TRUE)

## --- locate matrix-market triplet (names vary in the LCA zip) ---
ff <- list.files(IN_DIR, full.names = TRUE, recursive = TRUE)
mtxf <- ff[grepl("\\.mtx(\\.gz)?$", ff)]
barf <- ff[grepl("barcode.*\\.tsv(\\.gz)?$", basename(ff), ignore.case = TRUE)]
genf <- ff[grepl("(genes|features).*\\.tsv(\\.gz)?$", basename(ff), ignore.case = TRUE)]
stopifnot("no .mtx found in IN_DIR"     = length(mtxf) >= 1,
          "no barcodes .tsv found"      = length(barf) >= 1,
          "no genes/features .tsv found"= length(genf) >= 1)
mtxf <- mtxf[which.max(file.info(mtxf)$size)]
cat("Matrix :", basename(mtxf), "\nBarcodes:", basename(barf[1]), "\nFeatures:", basename(genf[1]), "\n")
M <- tryCatch(ReadMtx(mtx = mtxf, cells = barf[1], features = genf[1], feature.column = 1),
              error = function(e) ReadMtx(mtx = mtxf, cells = barf[1], features = genf[1], feature.column = 2))
cat("Loaded matrix:", nrow(M), "genes x", ncol(M), "cells\n")

## --- annotation ---
an <- as.data.frame(fread(ANNOT))
cat("\nAnnotation columns:", paste(colnames(an), collapse = ", "), "\n")
print(utils::head(an, 3))
## the cell-barcode column = the one whose values best match the matrix colnames
cell_col <- names(which.max(vapply(an, function(v) mean(as.character(v) %in% colnames(M)), numeric(1))))
cat("Cell-barcode column detected:", cell_col, "\n")
pick <- function(cands, override) if (!is.na(override)) override else {
 hit <- names(an)[tolower(names(an)) %in% cands]; if (length(hit)) hit[1] else NA }
donorc <- pick(c("patient","sample","subject","donor","patientid","sampleid"), DONOR_COL)
groupc <- pick(c("disease","condition","diet","group","typesample","status","diagnosis"), GROUP_COL)
stopifnot("could not find a donor column – set DONOR_COL manually" = !is.na(donorc))
cat("Using donor column:", donorc, "| group column:", if (is.na(groupc)) "(none -> 'all')" else groupc, "\n")

rownames(an) <- as.character(an[[cell_col]])
common <- intersect(colnames(M), rownames(an))
cat("Cells with annotation:", length(common), "/", ncol(M), "\n")
M <- M[, common]; an <- an[common, ]

obj <- CreateSeuratObject(M, project = "Guilliams2022", min.cells = 3, min.features = MIN_GENES)
obj$donor <- as.character(an[[donorc]])
obj$group <- if (is.na(groupc)) "all" else as.character(an[[groupc]])
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA >= MIN_GENES & nFeature_RNA <= MAX_GENES & percent.mt <= MAX_MT)

cat("\nAfter QC:", ncol(obj), "cells |", length(unique(obj$donor)), "donors\n")
cat("Donors by group:\n"); print(tapply(obj$donor, obj$group, function(x) length(unique(x))))
saveRDS(obj, OUT_RDS, compress = FALSE)
cat("\nWrote", OUT_RDS, "\n")
cat("Next: run sc_guilliams_clean.R (patient-level, diet-stratified Guilliams analysis).\n")
