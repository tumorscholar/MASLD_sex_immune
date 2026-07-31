## geo_loaders.R --------------------------------------------------------------
## Shared GEO download/parse helpers (Bioconductor GEOquery). Sourced by 01, 02 and 05.
## Uses Bioconductor's GEOquery. Metadata + array data are handled cleanly by
## GEOquery; RNA-seq COUNT tables live in each series' supplementary files, so
## for those we download the suppl files and apply simple format heuristics
## (largest count/tpm/matrix table; or per-sample merge when a
## series ships one file per GSM).
## ---------------------------------------------------------------------------
suppressMessages(library(GEOquery))
options(timeout = 1200)

## ---- per-sample characteristics -> named list ----
gsm_chars <- function(gsm) {
  ch <- Meta(gsm)[["characteristics_ch1"]]
  out <- list()
  if (is.null(ch)) return(out)
  for (c in ch) if (grepl(":", c)) {
    kv <- strsplit(c, ":", fixed = TRUE)[[1]]
    out[[tolower(trimws(kv[1]))]] <- trimws(paste(kv[-1], collapse = ":"))
  }
  out
}
## return the value of the first characteristic matching the HIGHEST-priority
## pattern (patterns are tried in order across all keys, so pattern order = priority)
find_val <- function(ch, patterns) {
  for (p in patterns) for (k in names(ch)) if (grepl(p, k)) return(ch[[k]])
  NA_character_
}
rec_sex <- function(ch) {
  v <- find_val(ch, c("\\bsex\\b", "gender"))
  if (is.na(v)) return(NA_character_)
  v <- tolower(trimws(v))
  if (startsWith(v, "m")) "M" else if (startsWith(v, "f")) "F" else NA_character_
}
rec_stage <- function(ch) {
  v <- find_val(ch, c("fibros")); if (is.na(v)) v <- find_val(ch, c("\\bstage\\b"))
  if (is.na(v)) {
    ## prefer a disease SUBTYPE (NAFL/NASH) over a constant disease-state field,
    ## so cohorts staged by subtype rather than fibrosis still get a severity axis
    g <- find_val(ch, c("subtype","group","diagnos","disease","condition"))
    return(if (!is.na(g)) paste0("grp:", g) else NA_character_)
  }
  m <- regmatches(v, regexpr("[Ff]?\\s*([0-4])", v))
  if (length(m) && nchar(m)) paste0("F", gsub("[^0-4]", "", m)) else trimws(v)
}
rec_num <- function(ch, patterns) {
  v <- find_val(ch, patterns); if (is.na(v)) return(NA_real_)
  m <- regmatches(v, regexpr("[-+]?\\d*\\.?\\d+", v))
  if (length(m) && nchar(m)) as.numeric(m) else NA_real_
}
rec_t2d <- function(ch) {
  v <- find_val(ch, c("diabet","t2dm?","\\bdm\\b")); if (is.na(v)) return(NA_real_)
  v <- tolower(trimws(v))
  if (v %in% c("1","yes","y","t2d","t2dm","diabetic","present","positive")) return(1)
  if (v %in% c("0","no","n","none","absent","negative","non-diabetic")) return(0)
  NA_real_
}
## flag non-baseline samples in longitudinal / interventional cohorts so they can
## be dropped (keep one cross-sectional observation per patient; e.g. GSE48452's
## post-bariatric-surgery follow-ups). Cross-sectional cohorts return FALSE for all.
rec_followup <- function(ch) {
  hay <- tolower(paste(unlist(ch), collapse = " | "))
  grepl("after surgery|post[- ]?surger|follow[- ]?up|post[- ]?treatment|post[- ]?intervention", hay)
}

## ---- array cohort: symbol x sample matrix from the ExpressionSet ----
## Collapses probes to gene symbols by mean.
array_symbol_matrix <- function(eset) {
  X <- exprs(eset)
  fd <- fData(eset)
  ## 1) prefer a clean gene-symbol column (Illumina, Affymetrix 3'-IVT like U133)
  sc <- grep("gene[ _]?symbol|^symbol$|ILMN_Gene", names(fd), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(sc)) {
    sym <- toupper(vapply(strsplit(as.character(fd[[sc]]), "[ /|]"), `[`, character(1), 1))
  } else {
    ## 2) Affymetrix Gene/Exon ST arrays (e.g. GPL11532, GPL16686): the symbol is
    ##    the 2nd '//' field of gene_assignment, e.g. "NM_x // SYMBOL // desc // ...".
    ga <- grep("gene[ _]?assignment", names(fd), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(ga)) stop("no gene-symbol / gene_assignment column in featureData; have: ",
                        paste(head(names(fd)), collapse = ", "))
    sym <- toupper(vapply(strsplit(as.character(fd[[ga]]), "\\s*//\\s*"),
                          function(p) if (length(p) >= 2) trimws(p[2]) else NA_character_,
                          character(1)))
  }
  keep <- !is.na(sym) & sym != "" & sym != "NA" & sym != "---"
  X <- X[keep, , drop = FALSE]; sym <- sym[keep]
  RS <- rowsum(X, group = sym)
  G  <- RS / as.vector(table(sym)[rownames(RS)])   # probe -> symbol mean
  G  # symbols x samples
}

## ---- rna-seq cohort: download suppl count tables, return native-id x sample ----
sep_sniff <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  head <- readChar(con, 8192, useBytes = TRUE); close(con)
  if (lengths(regmatches(head, gregexpr("\t", head))) >=
      lengths(regmatches(head, gregexpr(",",  head)))) "\t" else ","
}
read_table_any <- function(path) {
  s <- sep_sniff(path)
  data.table::fread(path, sep = s, header = TRUE, data.table = FALSE)
}
rnaseq_matrix <- function(gse_id) {
  ddir <- file.path(GEO_CACHE, gse_id)
  if (!dir.exists(ddir)) dir.create(ddir, recursive = TRUE, showWarnings = FALSE)
  tab_pat <- "\\.(csv|tsv|txt|tab)(\\.gz)?$"
  ## (re)download supplementary files whenever the cache has no tabular table yet
  ## -- an empty/partial folder from a failed prior run must NOT skip the download.
  if (!length(list.files(ddir, pattern = tab_pat, recursive = TRUE, ignore.case = TRUE)))
    getGEOSuppFiles(gse_id, baseDir = GEO_CACHE, makeDirectory = TRUE)
  ## untar any tarballs
  for (f in list.files(ddir, pattern = "\\.tar$", full.names = TRUE))
    try(untar(f, exdir = file.path(ddir, "untar")), silent = TRUE)
  files <- list.files(ddir, pattern = tab_pat,
                      recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("no tabular suppl for ", gse_id)
  if (length(files) > 20) {                       # per-sample layout (one file per GSM)
    ## if some files are GSM-named and some are not, the non-GSM ones are usually
    ## a bundled series-level matrix (e.g. GSE167523_Raw_gene_counts_matrix.txt.gz)
    ## that would otherwise be read in as a spurious extra "sample" -- drop them.
    has_gsm <- grepl("GSM\\d+", basename(files))
    if (any(has_gsm) && any(!has_gsm)) files <- files[has_gsm]
    cols <- list()
    for (f in files) {
      key <- regmatches(basename(f), regexpr("GSM\\d+", basename(f)))
      if (!length(key)) key <- basename(f)
      df <- try(read_table_any(f), silent = TRUE)
      if (!inherits(df, "try-error") && ncol(df) >= 2) {
        v <- suppressWarnings(as.numeric(df[[2]])); names(v) <- as.character(df[[1]])
        cols[[key]] <- v
      }
    }
    ids <- Reduce(union, lapply(cols, names))
    M <- sapply(cols, function(v) v[ids]); rownames(M) <- ids
    return(as.data.frame(M))
  }
  ## single matrix file: pick the largest count/tpm/matrix table
  cand <- files[grepl("count|fpkm|tpm|expr|matrix|norm|gene", basename(files), ignore.case = TRUE)]
  if (!length(cand)) cand <- files
  f <- cand[which.max(file.info(cand)$size)]
  df <- read_table_any(f); rn <- as.character(df[[1]]); df[[1]] <- NULL
  M <- as.data.frame(df); rownames(M) <- rn
  M
}

## map RNA-seq matrix column names (titles or GSMs) to GSM accessions
map_cols_to_gsm <- function(gse, M) {
  gsms <- GSMList(gse)
  ids <- names(gsms)
  t2g <- list()
  for (n in ids) for (t in Meta(gsms[[n]])[["title"]]) t2g[[trimws(t)]] <- n
  newc <- vapply(colnames(M), function(c) {
    if (c %in% ids) return(c)
    if (!is.null(t2g[[c]])) return(t2g[[c]])
    hit <- names(t2g)[vapply(names(t2g), function(t) grepl(c, t, fixed=TRUE) || grepl(t, c, fixed=TRUE), logical(1))]
    if (length(hit) == 1) t2g[[hit]] else c
  }, character(1))
  colnames(M) <- newc; M
}
