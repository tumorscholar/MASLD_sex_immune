## 02_build_matrix.R ----------------------------------------------------------
## Build the analysis matrix: one singscore per cell-type / functional-state
## signature per sample, across all cohorts, plus covariates. Genes are matched
## whether the cohort is indexed by symbol, Ensembl or Entrez.
## Writes analysis_matrix.csv (the input to 03/08/09) and signature_defs.json.
## ---------------------------------------------------------------------------
source("00_config.R")
source("geo_loaders.R")
suppressMessages({library(AnnotationDbi); library(org.Hs.eg.db); library(singscore)})

## Score a signature with the Bioconductor singscore package (Foroutan et al., 2018).
## Filters gene IDs to those present in the ranked matrix and returns the per-sample
## TotalScore; returns NULL if the up-set has no genes present (matches prior loop).
score_sig <- function(ranks, up_ids, dn_ids = NULL) {
  up_ids <- unique(up_ids[!is.na(up_ids) & up_ids %in% rownames(ranks)])
  if (!length(up_ids)) return(NULL)
  dn_ids <- if (!is.null(dn_ids)) unique(dn_ids[!is.na(dn_ids) & dn_ids %in% rownames(ranks)]) else NULL
  sc <- if (!is.null(dn_ids) && length(dn_ids))
          singscore::simpleScore(ranks, upSet = up_ids, downSet = dn_ids)
        else
          singscore::simpleScore(ranks, upSet = up_ids)
  sc$TotalScore
}

## ---- symbol -> Ensembl / Entrez maps for all signature genes ----
sym2ens <- tryCatch(mapIds(org.Hs.eg.db, ALL_SYMS, "ENSEMBL", "SYMBOL", multiVals = "first"),
                    error = function(e) setNames(rep(NA, length(ALL_SYMS)), ALL_SYMS))
sym2ent <- tryCatch(mapIds(org.Hs.eg.db, ALL_SYMS, "ENTREZID", "SYMBOL", multiVals = "first"),
                    error = function(e) setNames(rep(NA, length(ALL_SYMS)), ALL_SYMS))
cat(sprintf("  mapped %d/%d symbols to Ensembl, %d to Entrez\n",
            sum(!is.na(sym2ens)), length(ALL_SYMS), sum(!is.na(sym2ent))))

## build an uppercase + ensembl-stripped lookup of a matrix's row IDs
build_lookup <- function(M) {
  idx <- rownames(M); norm <- new.env(hash = TRUE)
  for (i in idx) {
    u <- toupper(i)
    if (is.null(norm[[u]])) assign(u, i, envir = norm)
    us <- sub("\\.\\d+$", "", u)
    if (is.null(norm[[us]])) assign(us, i, envir = norm)
  }
  norm
}
find_gene <- function(sym, norm) {
  for (cand in c(toupper(sym), toupper(sym2ens[[sym]] %||% ""), sym2ent[[sym]] %||% "")) {
    if (nzchar(cand) && !is.null(norm[[cand]])) return(norm[[cand]])
  }
  NA_character_
}
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

rows <- list()
for (co in COHORTS) {
  gid <- co$gse; typ <- co$type
  cat("\n======", gid, "======\n")
  res <- try({
    gse_soft <- getGEO_retry(gid, destdir = GEO_CACHE, GSEMatrix = FALSE)
    if (typ == "array") {
      eset <- getGEO_retry(gid, destdir = GEO_CACHE, GSEMatrix = TRUE, getGPL = TRUE)[[1]]
      M <- array_symbol_matrix(eset)
    } else {
      M <- as.matrix(map_cols_to_gsm(gse_soft, rnaseq_matrix(gid)))
    }
    M <- M[!duplicated(rownames(M)), , drop = FALSE]
    M <- M[rowSums(!is.na(M)) > 0, , drop = FALSE]
    M <- maybe_log2(M)
    norm <- build_lookup(M); ranks <- singscore::rankGenes(M)

    ## sex assignment (within-cohort z of XIST vs Y)
    gvec <- function(sym) { h <- find_gene(sym, norm); if (!is.na(h)) as.numeric(M[h, ]) else NULL }
    yv <- Filter(Negate(is.null), lapply(MALE_MARKERS, gvec))
    xv <- gvec("XIST")
    zY <- if (length(yv)) zscore(rowMeans(do.call(cbind, yv))) else rep(0, ncol(M))
    zX <- if (!is.null(xv)) zscore(xv) else rep(0, ncol(M))
    diff <- zY - zX
    sexs <- ifelse(diff > SEX_MARGIN, "M", ifelse(diff < -SEX_MARGIN, "F", "Ambiguous"))
    names(sexs) <- colnames(M)

    ## signature scores
    score_cols <- list()
    for (nm in names(CELLTYPE)) {
      ids <- vapply(CELLTYPE[[nm]], find_gene, character(1), norm = norm)
      s <- score_sig(ranks, ids); if (!is.null(s)) score_cols[[nm]] <- s
    }
    for (nm in names(STATE)) {
      up <- vapply(STATE[[nm]]$up, find_gene, character(1), norm = norm)
      dn <- if (length(STATE[[nm]]$down)) vapply(STATE[[nm]]$down, find_gene, character(1), norm = norm) else NULL
      s <- score_sig(ranks, up, dn); if (!is.null(s)) score_cols[[nm]] <- s
    }
    S <- as.data.frame(score_cols); rownames(S) <- colnames(M)

    ## covariates from the SOFT metadata
    gsms <- GSMList(gse_soft)
    for (n in rownames(S)) {
      if (!n %in% names(gsms)) next
      ch <- gsm_chars(gsms[[n]])
      if (isTRUE(rec_followup(ch))) next   # drop non-baseline (post-surgery/follow-up) samples
      r <- list(gsm = n, cohort = gid, platform = typ,
                sex_assigned = unname(sexs[n]), recorded_sex = rec_sex(ch),
                fibrosis_stage = rec_stage(ch), age = rec_num(ch, c("\\bage\\b")),
                bmi = rec_num(ch, c("\\bbmi\\b","body mass")), t2d = rec_t2d(ch))
      for (c in names(S)) r[[c]] <- S[n, c]
      rows[[length(rows)+1]] <- as.data.frame(r, stringsAsFactors = FALSE)
    }
    cat(sprintf("  scored %d samples | readouts: %d | sex F/M/Amb: %d/%d/%d\n",
        nrow(S), ncol(S), sum(sexs=="F"), sum(sexs=="M"), sum(sexs=="Ambiguous")))
  }, silent = TRUE)
  if (inherits(res, "try-error")) cat("  !! FAILED", gid, ":", conditionMessage(attr(res,"condition")), "\n")
}

df <- do.call(rbind, lapply(rows, function(r) {
  ## align columns across cohorts (some readouts may be missing in a cohort)
  r
}))
df <- as.data.frame(data.table::rbindlist(rows, fill = TRUE))
out <- file.path(REALDIR, "analysis_matrix.csv"); write.csv(df, out, row.names = FALSE)

## signature definitions (for provenance)
defs <- list(celltype = CELLTYPE,
             state = lapply(STATE, function(s) list(up = s$up, down = s$down)))
writeLines(jsonlite::toJSON(defs, pretty = TRUE, auto_unbox = TRUE),
           file.path(REALDIR, "signature_defs.json"))

cat("\nWrote", out, " shape", nrow(df), "x", ncol(df), "\n")
if (nrow(df)) {
  cat("Cohorts:", paste(names(table(df$cohort)), table(df$cohort), sep="=", collapse=" "), "\n")
  cat("Readout columns:", paste(grep("^(ct_|st_)", names(df), value = TRUE), collapse=", "), "\n")
}
