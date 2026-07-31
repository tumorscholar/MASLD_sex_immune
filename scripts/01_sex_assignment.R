## 01_sex_assignment.R --------------------------------------------------------
## Expression-based sex assignment for the public MASLD cohorts.
## Recovers sex from XIST (female) vs a 9-gene Y panel (male) by within-cohort
## z-scoring; diff = z_Y - z_XIST; call M if diff>0.5, F if diff< -0.5, else
## Ambiguous (dropped). Checks calls against recorded sex where available.
## ---------------------------------------------------------------------------
source("00_config.R")
source("geo_loaders.R")

OUTDIR <- file.path(REALDIR, "sex_assign_out"); dir.create(OUTDIR, showWarnings = FALSE)
ALL_MARK <- c(FEMALE_MARKERS, MALE_MARKERS)

## hardcoded cross-refs so we can find markers whether the matrix is indexed by
## symbol, Ensembl or Entrez.
ENSEMBL <- c(XIST="ENSG00000229807", RPS4Y1="ENSG00000129824", DDX3Y="ENSG00000067048",
  EIF1AY="ENSG00000198692", UTY="ENSG00000183878", KDM5D="ENSG00000012817",
  USP9Y="ENSG00000114374", NLGN4Y="ENSG00000165246", ZFY="ENSG00000067646",
  TXLNGY="ENSG00000131002")
ENTREZ <- c(XIST="7503", RPS4Y1="6192", DDX3Y="8653", EIF1AY="9086", UTY="7404",
  KDM5D="8284", USP9Y="8287", NLGN4Y="22829", ZFY="7544", TXLNGY="246126")

## pull marker rows out of a matrix whatever its ID type
gene_rows <- function(M) {
  M <- M[!duplicated(rownames(M)), , drop = FALSE]
  idx <- rownames(M)
  up  <- toupper(idx)
  ens <- toupper(sub("\\.\\d+$", "", idx))
  ent <- sub("\\.0$", "", idx)
  out <- list()
  for (g in ALL_MARK) {
    h <- NULL
    if (g %in% up) h <- idx[which(up == g)[1]]
    else if (ENSEMBL[[g]] %in% ens) h <- idx[which(ens == ENSEMBL[[g]])[1]]
    else if (ENTREZ[[g]] %in% ent) h <- idx[which(ent == ENTREZ[[g]])[1]]
    if (!is.null(h)) out[[g]] <- suppressWarnings(as.numeric(M[h, ]))
  }
  if (!length(out)) return(NULL)
  R <- as.data.frame(out); rownames(R) <- colnames(M)
  R  # samples x markers-present
}

classify_sex <- function(gr) {
  ym <- intersect(MALE_MARKERS, names(gr))
  hasX <- "XIST" %in% names(gr)
  if (!length(ym) && !hasX) stop("no XIST/Y markers present")
  zY <- if (length(ym)) zscore(rowMeans(gr[, ym, drop = FALSE], na.rm = TRUE)) else rep(0, nrow(gr))
  zX <- if (hasX) zscore(gr[["XIST"]]) else rep(0, nrow(gr))
  d  <- zY - zX
  call <- ifelse(d > SEX_MARGIN, "M", ifelse(d < -SEX_MARGIN, "F", "Ambiguous"))
  data.frame(gsm = rownames(gr), zXIST = zX, zYmean = zY, diff = d,
             sex_assigned = call, stringsAsFactors = FALSE)
}

allc <- list()
for (co in COHORTS) {
  gid <- co$gse; typ <- co$type
  cat("\n======", gid, sprintf("(%s)", typ), "======\n")
  res <- try({
    gse <- getGEO_retry(gid, destdir = GEO_CACHE, GSEMatrix = TRUE, getGPL = (typ == "array"))
    eset <- if (is.list(gse)) gse[[1]] else gse
    gsms <- GSMList(getGEO_retry(gid, destdir = GEO_CACHE, GSEMatrix = FALSE))
    meta <- do.call(rbind, lapply(names(gsms), function(n) {
      ch <- gsm_chars(gsms[[n]])
      data.frame(gsm = n, recorded_sex = rec_sex(ch),
                 fibrosis_stage = rec_stage(ch),
                 is_followup = rec_followup(ch), stringsAsFactors = FALSE)
    }))
    if (typ == "array") {
      G <- array_symbol_matrix(eset)                       # symbols x samples
      G <- maybe_log2(G); gr <- gene_rows(G)
    } else {
      M <- map_cols_to_gsm(getGEO_retry(gid, destdir = GEO_CACHE, GSEMatrix = FALSE),
                           rnaseq_matrix(gid))
      M <- as.matrix(M); M <- maybe_log2(M); gr <- gene_rows(M)
    }
    calls <- classify_sex(gr)
    df <- merge(meta, calls, by = "gsm"); df$gse <- gid
    if ("is_followup" %in% names(df)) df <- df[!df$is_followup, , drop = FALSE]  # baseline only
    kn <- df[!is.na(df$recorded_sex) & df$sex_assigned != "Ambiguous", ]
    if (nrow(kn))
      cat(sprintf("  Concordance vs recorded sex: %.1f%% (n=%d) | Ambiguous: %d\n",
          100*mean(kn$recorded_sex == kn$sex_assigned), nrow(kn),
          sum(df$sex_assigned == "Ambiguous")))
    else
      cat("  No recorded sex in GEO -> assignment only. Ambiguous:",
          sum(df$sex_assigned == "Ambiguous"), "\n")
    print(table(df$fibrosis_stage, df$sex_assigned))
    write.csv(table(df$fibrosis_stage, df$sex_assigned),
              file.path(OUTDIR, paste0("sex_by_stage_", gid, ".csv")))
    allc[[gid]] <- df
  }, silent = TRUE)
  if (inherits(res, "try-error")) cat("  !! FAILED", gid, ":", conditionMessage(attr(res,"condition")), "\n")
}

if (length(allc)) {
  A <- do.call(rbind, allc)
  write.csv(A, file.path(OUTDIR, "per_sample_sex_calls.csv"), row.names = FALSE)
  pooled <- table(paste(A$gse, A$fibrosis_stage), A$sex_assigned)
  write.csv(pooled, file.path(OUTDIR, "sex_by_stage_ALL.csv"))
  cat("\n===== POOLED sex x fibrosis-stage =====\n"); print(pooled)
  cat("\nWrote", OUTDIR, "-> sex_by_stage_ALL.csv\n")
} else cat("\nNo cohorts succeeded.\n")
