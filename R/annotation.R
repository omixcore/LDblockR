.parse_gff_attributes <- function(x) {
  parts <- strsplit(x, ";", fixed = TRUE)[[1L]]
  out <- character()
  for (part in parts) {
    part <- trimws(part)
    if (!nzchar(part)) next
    if (grepl("=", part, fixed = TRUE)) {
      z <- strsplit(part, "=", fixed = TRUE)[[1L]]
      key <- z[1L]
      value <- paste(z[-1L], collapse = "=")
    } else {
      z <- strsplit(part, "[[:space:]]+", perl = TRUE)[[1L]]
      key <- z[1L]
      value <- paste(z[-1L], collapse = " ")
    }
    value <- gsub('^"|"$', "", value)
    out[key] <- utils::URLdecode(value)
  }
  out
}

#' Read regional GFF3/GTF annotation
#'
#' @param file GFF3, GTF, or compressed annotation file.
#' @param region Optional genomic region.
#' @param feature_types Optional feature types to retain.
#' @return A normalized annotation data frame.
#' @export
read_gff3 <- function(file, region = NULL,
                      feature_types = c("gene", "mRNA", "transcript", "exon", "CDS",
                                        "five_prime_UTR", "three_prime_UTR", "UTR")) {
  reg <- parse_region(region)
  con <- .open_text(file)
  on.exit(close(con), add = TRUE)
  rows <- list()
  k <- 0L
  done <- FALSE
  while (!done) {
    lines <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    fields <- strsplit(lines, "\t", fixed = TRUE)
    for (f in fields) {
      if (length(f) < 9L) next
      chr <- f[1L]
      start <- suppressWarnings(as.numeric(f[4L]))
      end <- suppressWarnings(as.numeric(f[5L]))
      if (!is.finite(start) || !is.finite(end)) next
      if (!is.null(reg)) {
        if (chr != reg$chr[1L] || end < reg$start[1L]) next
        if (chr == reg$chr[1L] && start > reg$end[1L]) {
          done <- TRUE
          break
        }
      }
      if (!is.null(feature_types) && !f[3L] %in% feature_types) next
      a <- .parse_gff_attributes(f[9L])
      pick <- function(keys) {
        z <- a[keys]
        z <- z[!is.na(z) & nzchar(z)]
        if (length(z)) unname(z[1L]) else NA_character_
      }
      k <- k + 1L
      rows[[k]] <- data.frame(
        chr = chr, source = f[2L], type = f[3L], start = start, end = end,
        score = suppressWarnings(as.numeric(f[6L])), strand = f[7L], phase = f[8L],
        id = pick(c("ID", "transcript_id", "gene_id")),
        parent = pick(c("Parent", "gene_id")),
        name = pick(c("Name", "gene_name", "gene", "locus_tag", "ID")),
        attributes = f[9L], stringsAsFactors = FALSE
      )
    }
  }
  if (!k) {
    return(data.frame(chr = character(), source = character(), type = character(),
                      start = numeric(), end = numeric(), score = numeric(),
                      strand = character(), phase = character(), id = character(),
                      parent = character(), name = character(), attributes = character(),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  out <- out[order(.chromosome_rank(out$chr), out$chr, out$start, out$end), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("ld_genes", class(out))
  out
}

.has_header <- function(file) {
  con <- .open_text(file)
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) return(FALSE)
    if (nzchar(line) && !startsWith(trimws(line), "#")) break
  }
  z <- strsplit(trimws(line), "[[:space:],]+", perl = TRUE)[[1L]]
  any(grepl("^(chr|chrom|chromosome|pos|position|bp|p|pval|pvalue|logp)$", tolower(z)))
}

#' Read regional association statistics
#'
#' @param x File or data frame.
#' @param region Optional region.
#' @param chr_col,pos_col,p_col,id_col Explicit column names; auto-detected by default.
#' @param value_is_logp Whether the selected value is already -log10(P).
#' @param header Whether a file has a header; auto-detected by default.
#' @param sep File separator; empty means arbitrary whitespace.
#' @return A normalized data frame with `chr`, `pos`, `p`, `logp`, and `id`.
#' @details A bare `regional_gwas.tsv` filename resolves to the bundled
#' regional example after installation. Paths that include a directory are
#' treated as user-supplied paths.
#' @export
read_gwas <- function(x, region = NULL, chr_col = NULL, pos_col = NULL,
                      p_col = NULL, id_col = NULL, value_is_logp = FALSE,
                      header = NULL, sep = "") {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      # Resolve a bare filename against the package's installed example data.
      # Paths containing a directory are never redirected, so a typo in a
      # user-supplied path still produces the usual missing-file error.
      bundled <- system.file("extdata", basename(x), package = "LDblockR")
      if (identical(dirname(x), ".") && nzchar(bundled) && file.exists(bundled)) {
        x <- bundled
      } else {
        .stopf("GWAS file does not exist: %s. Use example_data(\"regional\") for the built-in regional example.", x)
      }
    }
    if (is.null(header)) header <- .has_header(x)
    con <- .open_text(x)
    on.exit(close(con), add = TRUE)
    tab <- utils::read.table(con, header = header, sep = sep, quote = "",
                             comment.char = "#", check.names = FALSE,
                             stringsAsFactors = FALSE)
    if (!header) {
      if (ncol(tab) < 3L) .stopf("Headerless GWAS input needs at least chr, position, and P/value columns.")
      names(tab)[1:3] <- c("chr", "pos", "p")
    }
  } else {
    tab <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  if (!nrow(tab)) .stopf("GWAS table is empty.")
  chr_col <- .null_coalesce(chr_col, .match_column(tab, c("chr", "chrom", "chromosome", "#chrom"), TRUE, "chromosome"))
  pos_col <- .null_coalesce(pos_col, .match_column(tab, c("pos", "position", "bp", "site"), TRUE, "position"))
  if (is.null(p_col)) {
    p_col <- .match_column(tab, c("p", "pval", "pvalue", "p_value", "p.value", "logp", "neglog10p", "value"), TRUE, "P/value")
    if (tolower(p_col) %in% c("logp", "neglog10p")) value_is_logp <- TRUE
  }
  pve_col <- .match_column(tab, c("PVE", "pve", "var_explained", "variance_explained",
                                 "percent_variance_explained"), FALSE)
  beta_col <- .match_column(tab, c("beta", "effect", "estimate", "coefficient"), FALSE)
  se_col <- .match_column(tab, c("se", "stderr", "standard_error"), FALSE)
  if (is.null(id_col)) id_col <- .match_column(tab, c("id", "snp", "marker", "rsid", "name"), FALSE)
  value <- as.numeric(tab[[p_col]])
  if (value_is_logp) {
    logp <- value
    p <- 10^(-value)
  } else {
    p <- value
    logp <- .safe_log10(p)
  }
  out <- data.frame(
    chr = as.character(tab[[chr_col]]), pos = as.numeric(tab[[pos_col]]),
    p = p, logp = logp,
    id = if (is.null(id_col)) paste0(tab[[chr_col]], ":", tab[[pos_col]]) else as.character(tab[[id_col]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(pve_col)) out$PVE <- as.numeric(tab[[pve_col]])
  if (!is.null(beta_col)) out$beta <- as.numeric(tab[[beta_col]])
  if (!is.null(se_col)) out$se <- as.numeric(tab[[se_col]])
  out <- out[is.finite(out$pos) & is.finite(out$logp) & .in_region(out$chr, out$pos, parse_region(region)), , drop = FALSE]
  out <- out[order(.chromosome_rank(out$chr), out$chr, out$pos), , drop = FALSE]
  rownames(out) <- NULL
  out
}
