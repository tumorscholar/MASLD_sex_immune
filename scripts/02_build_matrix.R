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
cov_rows <- list()   # per-cohort x per-signature marker-gene coverage (QC)
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
    norm <- build_lookup(M); ranks <- singscore::rankGenes(M, tiesMethod = "average")
    
    ## sex assignment (within-cohort z of XIST vs Y)
    gvec <- function(sym) { h <- find_gene(sym, norm); if (!is.na(h)) as.numeric(M[h, ]) else NULL }
    yv <- Filter(Negate(is.null), lapply(MALE_MARKERS, gvec))
    xv <- gvec("XIST")
    zY <- if (length(yv)) zscore(rowMeans(do.call(cbind, yv))) else rep(0, ncol(M))
    zX <- if (!is.null(xv)) zscore(xv) else rep(0, ncol(M))
    diff <- zY - zX
    sexs <- ifelse(diff > SEX_MARGIN, "M", ifelse(diff < -SEX_MARGIN, "F", "Ambiguous"))
    names(sexs) <- colnames(M)
    
    ## record marker-gene coverage of each signature in this cohort (QC)
    log_cov <- function(nm, genes, ids) {
      miss <- genes[is.na(ids)]
      cov_rows[[length(cov_rows) + 1]] <<- data.frame(
        cohort = gid, platform = typ, signature = nm,
        n_total = length(genes), n_detected = sum(!is.na(ids)),
        coverage = round(sum(!is.na(ids)) / length(genes), 3),
        missing = paste(miss, collapse = ";"), stringsAsFactors = FALSE)
    }
    
    ## signature scores
    score_cols <- list()
    for (nm in names(CELLTYPE)) {
      ids <- vapply(CELLTYPE[[nm]], find_gene, character(1), norm = norm)
      log_cov(nm, CELLTYPE[[nm]], ids)
      s <- score_sig(ranks, ids); if (!is.null(s)) score_cols[[nm]] <- s
    }
    for (nm in names(STATE)) {
      up <- vapply(STATE[[nm]]$up, find_gene, character(1), norm = norm)
      dn <- if (length(STATE[[nm]]$down)) vapply(STATE[[nm]]$down, find_gene, character(1), norm = norm) else NULL
      allg <- c(STATE[[nm]]$up, if (length(STATE[[nm]]$down)) STATE[[nm]]$down else character(0))
      allids <- c(up, if (!is.null(dn)) dn else character(0))
      log_cov(nm, allg, allids)
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

## ---- signature coverage QC ------------------------------------------------
## How many marker genes of each signature were detected in each cohort. Long
## form (one row per cohort x signature) plus a compact detected/total matrix.
if (length(cov_rows)) {
  cov <- as.data.frame(data.table::rbindlist(cov_rows, fill = TRUE))
  write.csv(cov, file.path(REALDIR, "signature_coverage.csv"), row.names = FALSE)
  cov$cell <- paste0(cov$n_detected, "/", cov$n_total)
  wide <- reshape(cov[, c("signature", "cohort", "cell")],
                  idvar = "signature", timevar = "cohort", direction = "wide")
  names(wide) <- sub("^cell\\.", "", names(wide))
  write.csv(wide, file.path(REALDIR, "signature_coverage_matrix.csv"), row.names = FALSE)
  low <- cov[cov$coverage < 0.5, ]
  cat(sprintf("\nSignature coverage: %d cohort x signature cells | median coverage %.0f%% | %d below 50%%\n",
              nrow(cov), 100 * median(cov$coverage), nrow(low)))
  if (nrow(low)) {
    cat("  Below 50% coverage (interpret these readouts cautiously):\n")
    for (i in seq_len(nrow(low)))
      cat(sprintf("    %-14s %-16s %d/%d  missing: %s\n",
                  low$cohort[i], low$signature[i], low$n_detected[i], low$n_total[i], low$missing[i]))
  }
  cat("Wrote signature_coverage.csv and signature_coverage_matrix.csv\n")
}

cat("\nWrote", out, " shape", nrow(df), "x", ncol(df), "\n")
if (nrow(df)) {
  cat("Cohorts:", paste(names(table(df$cohort)), table(df$cohort), sep="=", collapse=" "), "\n")
  cat("Readout columns:", paste(grep("^(ct_|st_)", names(df), value = TRUE), collapse=", "), "\n")
}
