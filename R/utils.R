.stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
.warnf <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)

.null_coalesce <- function(x, y) if (is.null(x)) y else x

.match_column <- function(x, candidates, required = TRUE, label = NULL) {
  nms <- names(x)
  hit <- match(tolower(candidates), tolower(nms), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) return(nms[hit[1L]])
  if (required) {
    .stopf("Cannot find %s column. Accepted names: %s",
           .null_coalesce(label, "required"), paste(candidates, collapse = ", "))
  }
  NULL
}

.read_sample_ids <- function(samples) {
  if (is.null(samples)) return(NULL)
  if (length(samples) == 1L && is.character(samples) && file.exists(samples)) {
    samples <- scan(samples, what = character(), quiet = TRUE)
  }
  unique(as.character(samples[nzchar(as.character(samples))]))
}

.open_text <- function(file) {
  if (!file.exists(file)) .stopf("File does not exist: %s", file)
  if (grepl("\\.(gz|bgz)$", file, ignore.case = TRUE)) {
    gzfile(file, open = "rt")
  } else if (grepl("\\.bz2$", file, ignore.case = TRUE)) {
    bzfile(file, open = "rt")
  } else if (grepl("\\.xz$", file, ignore.case = TRUE)) {
    xzfile(file, open = "rt")
  } else {
    base::file(file, open = "rt")
  }
}

#' Parse a genomic region
#'
#' @param region A string such as `chr1:1000-2000` or `chr1:1000:2000`.
#' @return A one-row data frame with `chr`, `start`, and `end`.
#' @export
parse_region <- function(region) {
  if (is.null(region) || !length(region) || is.na(region[1L]) || !nzchar(region[1L])) {
    return(NULL)
  }
  region <- gsub(",", "", trimws(as.character(region[1L])), fixed = TRUE)
  z <- regexec("^([^:]+):([0-9]+)(?:-|:)([0-9]+)$", region, perl = TRUE)
  p <- regmatches(region, z)[[1L]]
  if (length(p) != 4L) {
    .stopf("Invalid region '%s'. Use chr:start-end (for example chr1:1000-2000).", region)
  }
  start <- as.numeric(p[3L])
  end <- as.numeric(p[4L])
  if (!is.finite(start) || !is.finite(end) || start < 0 || end < start) {
    .stopf("Invalid region coordinates in '%s'.", region)
  }
  data.frame(chr = p[2L], start = start, end = end, stringsAsFactors = FALSE)
}

.chr_key <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z[is.na(x)] <- NA_character_
  sub("^chr", "", z)
}

.same_chr <- function(x, y) {
  lhs <- .chr_key(x)
  rhs <- .chr_key(y)
  !is.na(lhs) & !is.na(rhs) & lhs == rhs
}

.in_region <- function(chr, pos, region) {
  if (is.null(region)) return(rep(TRUE, length(pos)))
  .same_chr(chr, region$chr[1L]) & pos >= region$start[1L] & pos <= region$end[1L]
}

.chromosome_rank <- function(chr) {
  x <- sub("^chr", "", as.character(chr), ignore.case = TRUE)
  num <- suppressWarnings(as.numeric(x))
  known <- c(X = 1e6 + 1, Y = 1e6 + 2, Z = 1e6 + 3, W = 1e6 + 4,
             MT = 1e6 + 5, M = 1e6 + 5)
  out <- num
  idx <- is.na(out) & toupper(x) %in% names(known)
  out[idx] <- known[toupper(x[idx])]
  idx <- is.na(out)
  if (any(idx)) out[idx] <- 2e6 + match(x[idx], sort(unique(x[idx])))
  out
}

.safe_log10 <- function(x) {
  x <- as.numeric(x)
  finite_positive <- x[is.finite(x) & x > 0]
  floor_value <- if (length(finite_positive)) min(finite_positive) / 10 else .Machine$double.xmin
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  out[ok] <- -log10(pmax(x[ok], floor_value))
  out
}

.hwe_p <- function(g) {
  g <- g[!is.na(g)]
  if (length(g) < 3L) return(NA_real_)
  obs <- tabulate(as.integer(g) + 1L, nbins = 3L)
  n <- sum(obs)
  p <- (obs[2L] + 2 * obs[3L]) / (2 * n)
  expct <- c(n * (1 - p)^2, n * 2 * p * (1 - p), n * p^2)
  if (any(expct <= 0)) return(NA_real_)
  stats::pchisq(sum((obs - expct)^2 / expct), df = 1, lower.tail = FALSE)
}

.variant_statistics <- function(genotypes) {
  n <- nrow(genotypes)
  called <- colSums(!is.na(genotypes))
  alt_count <- colSums(genotypes, na.rm = TRUE)
  af <- alt_count / (2 * called)
  af[called == 0L] <- NA_real_
  data.frame(
    af = af,
    maf = pmin(af, 1 - af),
    missing_rate = 1 - called / n,
    het_rate = colSums(genotypes == 1, na.rm = TRUE) / pmax(called, 1L),
    call_count = called,
    hwe_p = vapply(seq_len(ncol(genotypes)), function(j) .hwe_p(genotypes[, j]), numeric(1L)),
    stringsAsFactors = FALSE
  )
}

.empty_blocks <- function() {
  out <- data.frame(
    block = character(), chr = character(), start = numeric(), end = numeric(),
    length_bp = numeric(), n_snps = integer(), start_index = integer(),
    end_index = integer(), snps = character(), method = character(),
    strong_fraction = numeric(), stringsAsFactors = FALSE
  )
  class(out) <- c("ld_blocks", class(out))
  out
}

.normalize_map <- function(map, n_variants, ids = NULL) {
  if (is.null(map)) {
    ids <- .null_coalesce(ids, paste0("SNP", seq_len(n_variants)))
    return(data.frame(chr = "1", pos = seq_len(n_variants), id = ids,
                      ref = NA_character_, alt = NA_character_, stringsAsFactors = FALSE))
  }
  if (is.atomic(map) && !is.data.frame(map)) {
    if (length(map) != n_variants) .stopf("map has %d positions but genotype matrix has %d variants.", length(map), n_variants)
    ids <- .null_coalesce(ids, paste0("SNP", seq_len(n_variants)))
    return(data.frame(chr = "1", pos = as.numeric(map), id = ids,
                      ref = NA_character_, alt = NA_character_, stringsAsFactors = FALSE))
  }
  map <- as.data.frame(map, stringsAsFactors = FALSE)
  if (nrow(map) != n_variants) .stopf("map has %d rows but genotype matrix has %d variants.", nrow(map), n_variants)
  chr_col <- .match_column(map, c("chr", "chrom", "chromosome", "#chrom"), FALSE)
  pos_col <- .match_column(map, c("pos", "position", "bp", "site"), TRUE, "position")
  id_col <- .match_column(map, c("id", "snp", "marker", "rsid", "name"), FALSE)
  ref_col <- .match_column(map, c("ref", "a1", "allele1", "reference"), FALSE)
  alt_col <- .match_column(map, c("alt", "a2", "allele2", "alternate"), FALSE)
  out <- data.frame(
    chr = if (is.null(chr_col)) "1" else as.character(map[[chr_col]]),
    pos = as.numeric(map[[pos_col]]),
    id = if (is.null(id_col)) .null_coalesce(ids, paste0("SNP", seq_len(n_variants))) else as.character(map[[id_col]]),
    ref = if (is.null(ref_col)) NA_character_ else as.character(map[[ref_col]]),
    alt = if (is.null(alt_col)) NA_character_ else as.character(map[[alt_col]]),
    stringsAsFactors = FALSE
  )
  extra <- setdiff(names(map), c(chr_col, pos_col, id_col, ref_col, alt_col))
  if (length(extra)) out[extra] <- map[extra]
  out
}
