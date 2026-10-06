#' Construct an LD genotype object
#'
#' @param genotypes Numeric matrix coded 0/1/2, with samples in rows.
#' @param map Variant map containing at least a position column.
#' @param sample_ids Optional sample identifiers.
#' @param haplotypes Optional 0/1 matrix with two haplotype rows per sample.
#' @param source Optional source description.
#' @return An object of class `ld_data`.
#' @export
as_ld_data <- function(genotypes, map = NULL, sample_ids = NULL,
                       haplotypes = NULL, source = NULL) {
  genotypes <- as.matrix(genotypes)
  storage.mode(genotypes) <- "numeric"
  if (!length(genotypes) || nrow(genotypes) < 1L || ncol(genotypes) < 1L) {
    .stopf("genotypes must be a non-empty samples-by-variants matrix.")
  }
  bad <- !is.na(genotypes) & !(genotypes %in% c(0, 1, 2))
  if (any(bad)) .stopf("Genotypes must be hard calls coded 0, 1, 2, or NA.")
  if (is.null(sample_ids)) sample_ids <- rownames(genotypes)
  if (is.null(sample_ids)) sample_ids <- paste0("sample", seq_len(nrow(genotypes)))
  if (length(sample_ids) != nrow(genotypes)) .stopf("sample_ids length does not match genotype rows.")
  sample_ids <- make.unique(as.character(sample_ids))
  rownames(genotypes) <- sample_ids
  variants <- .normalize_map(map, ncol(genotypes), colnames(genotypes))
  if (any(!is.finite(variants$pos))) .stopf("Variant positions must be finite numbers.")
  variants$id[is.na(variants$id) | !nzchar(variants$id)] <- paste0(
    variants$chr[is.na(variants$id) | !nzchar(variants$id)], ":",
    variants$pos[is.na(variants$id) | !nzchar(variants$id)]
  )
  variants$id <- make.unique(variants$id)
  ord <- order(.chromosome_rank(variants$chr), variants$chr, variants$pos, variants$id)
  variants <- variants[ord, , drop = FALSE]
  genotypes <- genotypes[, ord, drop = FALSE]
  colnames(genotypes) <- variants$id
  if (!is.null(haplotypes)) {
    haplotypes <- as.matrix(haplotypes)
    storage.mode(haplotypes) <- "numeric"
    if (nrow(haplotypes) != 2L * nrow(genotypes) || ncol(haplotypes) != ncol(genotypes)) {
      .stopf("haplotypes must have two rows per sample and one column per variant.")
    }
    if (any(!is.na(haplotypes) & !(haplotypes %in% c(0, 1)))) {
      .stopf("haplotypes must be coded 0, 1, or NA.")
    }
    haplotypes <- haplotypes[, ord, drop = FALSE]
    colnames(haplotypes) <- variants$id
    rownames(haplotypes) <- as.vector(rbind(paste0(sample_ids, ".1"), paste0(sample_ids, ".2")))
  }
  stats <- .variant_statistics(genotypes)
  duplicate_stats <- intersect(names(variants), names(stats))
  if (length(duplicate_stats)) variants[duplicate_stats] <- NULL
  variants <- cbind(variants, stats)
  rownames(variants) <- NULL
  out <- list(
    genotypes = genotypes,
    variants = variants,
    samples = sample_ids,
    haplotypes = haplotypes,
    source = source,
    n_samples = nrow(genotypes),
    n_variants = ncol(genotypes)
  )
  class(out) <- "ld_data"
  out
}

#' @export
print.ld_data <- function(x, ...) {
  cat("<ld_data>", x$n_samples, "samples x", x$n_variants, "variants\n")
  cat("  region:", paste0(x$variants$chr[1L], ":", format(min(x$variants$pos), scientific = FALSE),
                           "-", format(max(x$variants$pos), scientific = FALSE)), "\n")
  cat("  mean MAF:", format(mean(x$variants$maf, na.rm = TRUE), digits = 3),
      " | mean missing:", format(mean(x$variants$missing_rate, na.rm = TRUE), digits = 3), "\n")
  if (!is.null(x$source)) cat("  source:", x$source, "\n")
  invisible(x)
}

#' Filter variants by standard quality statistics
#'
#' @param x An `ld_data` object.
#' @param min_maf Minimum minor allele frequency.
#' @param max_maf Maximum minor allele frequency.
#' @param max_missing Maximum missing-call fraction.
#' @param max_het Maximum heterozygous-call fraction.
#' @param min_hwe_p Minimum chi-square Hardy-Weinberg P value; use zero to disable.
#' @param quiet Suppress the filtering summary.
#' @return A filtered `ld_data` object.
#' @export
filter_variants <- function(x, min_maf = 0.05, max_maf = 0.5,
                            max_missing = 0.25, max_het = 1,
                            min_hwe_p = 0, quiet = FALSE) {
  if (!inherits(x, "ld_data")) .stopf("x must be an ld_data object.")
  v <- x$variants
  keep <- is.finite(v$maf) & v$maf >= min_maf & v$maf <= max_maf &
    v$missing_rate <= max_missing & v$het_rate <= max_het
  if (min_hwe_p > 0) keep <- keep & !is.na(v$hwe_p) & v$hwe_p >= min_hwe_p
  if (!any(keep)) {
    .stopf("No variants remain after filtering (MAF >= %g, missing <= %g, heterozygosity <= %g).",
           min_maf, max_missing, max_het)
  }
  out <- subset_ld_data(x, variants = keep)
  if (!quiet) message(sprintf("Retained %d of %d variants after filtering.", out$n_variants, x$n_variants))
  out
}

#' Subset an LD genotype object
#'
#' @param x An `ld_data` object.
#' @param variants Variant indices, IDs, or logical vector.
#' @param samples Sample indices, IDs, or logical vector.
#' @return A new `ld_data` object with recalculated variant statistics.
#' @export
subset_ld_data <- function(x, variants = NULL, samples = NULL) {
  if (!inherits(x, "ld_data")) .stopf("x must be an ld_data object.")
  resolve <- function(sel, ids, what) {
    if (is.null(sel)) return(seq_along(ids))
    if (is.character(sel)) {
      idx <- match(sel, ids)
      if (anyNA(idx)) .stopf("Unknown %s: %s", what, paste(sel[is.na(idx)], collapse = ", "))
      return(idx)
    }
    if (is.logical(sel)) {
      if (length(sel) != length(ids)) .stopf("Logical %s selector has the wrong length.", what)
      return(which(sel))
    }
    idx <- as.integer(sel)
    if (anyNA(idx) || any(idx < 1L | idx > length(ids))) .stopf("Invalid %s indices.", what)
    idx
  }
  vi <- resolve(variants, x$variants$id, "variant")
  si <- resolve(samples, x$samples, "sample")
  if (!length(vi) || !length(si)) .stopf("Subsetting cannot produce an empty genotype object.")
  hap <- NULL
  if (!is.null(x$haplotypes)) {
    hi <- as.vector(rbind(2L * si - 1L, 2L * si))
    hap <- x$haplotypes[hi, vi, drop = FALSE]
  }
  map_cols <- setdiff(names(x$variants), c("af", "maf", "missing_rate", "het_rate", "call_count", "hwe_p"))
  as_ld_data(x$genotypes[si, vi, drop = FALSE], x$variants[vi, map_cols, drop = FALSE],
             sample_ids = x$samples[si], haplotypes = hap, source = x$source)
}
