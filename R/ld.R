.pair_haplotype_stats <- function(g1 = NULL, g2 = NULL, h1 = NULL, h2 = NULL,
                                  tol = 1e-10, max_iter = 1000L) {
  if (!is.null(h1) && !is.null(h2)) {
    ok <- !is.na(h1) & !is.na(h2)
    h1 <- h1[ok]
    h2 <- h2[ok]
    if (length(h1) < 4L) return(list(n = floor(length(h1) / 2), freq = rep(NA_real_, 4L),
                                     pA = NA_real_, pB = NA_real_, D = NA_real_,
                                     dprime = NA_real_, r2 = NA_real_))
    counts <- tabulate(2L * as.integer(h1) + as.integer(h2) + 1L, nbins = 4L)
    freq <- counts / sum(counts)
    n <- floor(sum(counts) / 2)
  } else {
    ok <- !is.na(g1) & !is.na(g2)
    g1 <- as.integer(g1[ok])
    g2 <- as.integer(g2[ok])
    n <- length(g1)
    if (n < 2L) return(list(n = n, freq = rep(NA_real_, 4L), pA = NA_real_, pB = NA_real_,
                            D = NA_real_, dprime = NA_real_, r2 = NA_real_))
    idx <- 3L * g1 + g2 + 1L
    nm <- matrix(tabulate(idx, nbins = 9L), nrow = 3L, byrow = TRUE,
                 dimnames = list(0:2, 0:2))
    base <- c(
      2 * nm[1, 1] + nm[1, 2] + nm[2, 1],
      2 * nm[1, 3] + nm[1, 2] + nm[2, 3],
      2 * nm[3, 1] + nm[2, 1] + nm[3, 2],
      2 * nm[3, 3] + nm[2, 3] + nm[3, 2]
    )
    double_het <- nm[2, 2]
    pA <- mean(g1) / 2
    pB <- mean(g2) / 2
    freq <- c((1 - pA) * (1 - pB), (1 - pA) * pB, pA * (1 - pB), pA * pB)
    for (iter in seq_len(max_iter)) {
      denom <- freq[1L] * freq[4L] + freq[2L] * freq[3L]
      coupling <- if (denom <= 0) 0.5 else freq[1L] * freq[4L] / denom
      counts <- base + double_het * c(coupling, 1 - coupling, 1 - coupling, coupling)
      next_freq <- counts / (2 * n)
      if (max(abs(next_freq - freq)) < tol) {
        freq <- next_freq
        break
      }
      freq <- next_freq
    }
  }
  names(freq) <- c("00", "01", "10", "11")
  pA <- freq[3L] + freq[4L]
  pB <- freq[2L] + freq[4L]
  D <- freq[4L] - pA * pB
  if (D >= 0) {
    dmax <- min(pA * (1 - pB), (1 - pA) * pB)
  } else {
    dmax <- min(pA * pB, (1 - pA) * (1 - pB))
  }
  denom_r2 <- pA * (1 - pA) * pB * (1 - pB)
  list(
    n = n, freq = freq, pA = pA, pB = pB, D = D,
    dprime = if (dmax > 0) min(1, abs(D) / dmax) else NA_real_,
    r2 = if (denom_r2 > 0) min(1, D^2 / denom_r2) else NA_real_
  )
}

.genotype_pair_counts <- function(g1, g2) {
  ok <- !is.na(g1) & !is.na(g2)
  g1 <- as.integer(g1[ok])
  g2 <- as.integer(g2[ok])
  matrix(tabulate(3L * g1 + g2 + 1L, nbins = 9L), nrow = 3L, byrow = TRUE)
}

.dprime_profile_ci <- function(g1, g2, step = 0.01) {
  s <- .pair_haplotype_stats(g1, g2)
  if (!is.finite(s$dprime) || s$n < 3L) return(c(lower = NA_real_, upper = NA_real_))
  pA <- s$pA
  pB <- s$pB
  sign_d <- if (s$D < 0) -1 else 1
  dmax <- if (sign_d > 0) min(pA * (1 - pB), (1 - pA) * pB) else min(pA * pB, (1 - pA) * (1 - pB))
  grid <- seq(0, 1, by = step)
  if (tail(grid, 1L) < 1) grid <- c(grid, 1)
  nm <- .genotype_pair_counts(g1, g2)
  ll <- rep(-Inf, length(grid))
  for (k in seq_along(grid)) {
    D <- sign_d * grid[k] * dmax
    f00 <- (1 - pA) * (1 - pB) + D
    f01 <- (1 - pA) * pB - D
    f10 <- pA * (1 - pB) - D
    f11 <- pA * pB + D
    f <- pmax(c(f00, f01, f10, f11), 0)
    prob <- matrix(c(
      f[1]^2, 2 * f[1] * f[2], f[2]^2,
      2 * f[1] * f[3], 2 * (f[1] * f[4] + f[2] * f[3]), 2 * f[2] * f[4],
      f[3]^2, 2 * f[3] * f[4], f[4]^2
    ), nrow = 3L, byrow = TRUE)
    used <- nm > 0
    if (all(prob[used] > 0)) ll[k] <- sum(nm[used] * log(prob[used]))
  }
  if (!any(is.finite(ll))) return(c(lower = NA_real_, upper = NA_real_))
  w <- exp(ll - max(ll, na.rm = TRUE))
  w[!is.finite(w)] <- 0
  if (sum(w) <= 0) return(c(lower = NA_real_, upper = NA_real_))
  cw <- cumsum(w / sum(w))
  c(lower = grid[which(cw >= 0.05)[1L]], upper = grid[which(cw >= 0.95)[1L]])
}

.pair_distance_mask <- function(variants, max_distance) {
  same_chr <- outer(as.character(variants$chr), as.character(variants$chr), `==`)
  dist <- abs(outer(variants$pos, variants$pos, `-`))
  same_chr & dist <= max_distance
}

#' Calculate pairwise linkage disequilibrium
#'
#' @param x An `ld_data` object.
#' @param measure `r2`, `dprime`, or `both`.
#' @param r2_method Automatic selection, dosage correlation, or phased haplotype correlation.
#' @param max_distance Maximum physical distance in base pairs.
#' @param min_n Minimum pairwise sample count.
#' @param ci Calculate profile-likelihood 90 percent confidence intervals for D-prime.
#' @param ci_step D-prime likelihood grid step.
#' @param use_phased Use VCF phase information when complete for a variant pair.
#' @return An object of class `ld_result`.
#' @export
ld_compute <- function(x, measure = c("r2", "dprime", "both"),
                       r2_method = c("auto", "dosage", "haplotype"),
                       max_distance = Inf, min_n = 5L,
                       ci = FALSE, ci_step = 0.01, use_phased = TRUE) {
  if (!inherits(x, "ld_data")) .stopf("x must be an ld_data object.")
  measure <- match.arg(measure)
  r2_method <- match.arg(r2_method)
  if (x$n_variants < 2L) .stopf("At least two variants are required for LD calculation.")
  if (!is.finite(max_distance)) max_distance <- Inf
  mask <- .pair_distance_mask(x$variants, max_distance)
  called <- !is.na(x$genotypes)
  nmat <- crossprod(called)
  dimnames(nmat) <- list(x$variants$id, x$variants$id)
  need_r2 <- measure %in% c("r2", "both")
  need_d <- measure %in% c("dprime", "both") || ci
  r2 <- NULL
  if (need_r2) {
    phase_all <- !is.null(x$haplotypes) &&
      (!"phase_complete" %in% names(x$variants) || all(x$variants$phase_complete))
    if (r2_method == "auto") r2_method <- if (phase_all) "haplotype" else "dosage"
    source_matrix <- x$genotypes
    if (r2_method == "haplotype") {
      if (is.null(x$haplotypes)) .stopf("Phased haplotypes are unavailable; use r2_method='dosage'.")
      if (!phase_all) .stopf("Not all heterozygous genotypes are phased; use r2_method='dosage'.")
      source_matrix <- x$haplotypes
    }
    r2 <- suppressWarnings(stats::cor(source_matrix, use = "pairwise.complete.obs")^2)
    if (length(r2) == 1L) r2 <- matrix(r2, 1L, 1L)
    dimnames(r2) <- list(x$variants$id, x$variants$id)
    r2[!mask | nmat < min_n] <- NA_real_
    diag(r2) <- 1
  }
  dprime <- dcoef <- lower <- upper <- NULL
  if (need_d) {
    m <- x$n_variants
    dprime <- matrix(NA_real_, m, m, dimnames = list(x$variants$id, x$variants$id))
    dcoef <- dprime
    if (ci) {
      lower <- dprime
      upper <- dprime
    }
    diag(dprime) <- 1
    diag(dcoef) <- 0
    if (ci) {
      diag(lower) <- 1
      diag(upper) <- 1
    }
    phase_complete <- if ("phase_complete" %in% names(x$variants)) as.logical(x$variants$phase_complete) else rep(FALSE, m)
    for (i in seq_len(m - 1L)) {
      for (j in (i + 1L):m) {
        if (!mask[i, j] || nmat[i, j] < min_n) next
        use_h <- use_phased && !is.null(x$haplotypes) && isTRUE(phase_complete[i]) && isTRUE(phase_complete[j])
        s <- if (use_h) {
          .pair_haplotype_stats(h1 = x$haplotypes[, i], h2 = x$haplotypes[, j])
        } else {
          .pair_haplotype_stats(x$genotypes[, i], x$genotypes[, j])
        }
        dprime[i, j] <- dprime[j, i] <- s$dprime
        dcoef[i, j] <- dcoef[j, i] <- s$D
        if (ci) {
          z <- .dprime_profile_ci(x$genotypes[, i], x$genotypes[, j], ci_step)
          lower[i, j] <- lower[j, i] <- z[1L]
          upper[i, j] <- upper[j, i] <- z[2L]
        }
      }
    }
  }
  out <- list(data = x, r2 = r2, dprime = dprime, D = dcoef,
              n = nmat, ci_lower = lower, ci_upper = upper,
              measure = measure, r2_method = r2_method,
              max_distance = max_distance, min_n = min_n,
              call = match.call())
  class(out) <- "ld_result"
  out
}

#' @export
print.ld_result <- function(x, ...) {
  cat("<ld_result>", x$data$n_samples, "samples x", x$data$n_variants, "variants\n")
  available <- c(if (!is.null(x$r2)) "r2", if (!is.null(x$dprime)) "Dprime",
                 if (!is.null(x$ci_lower)) "Dprime-CI")
  cat("  statistics:", paste(available, collapse = ", "), "\n")
  cat("  maximum distance:", if (is.finite(x$max_distance)) paste(x$max_distance, "bp") else "unlimited", "\n")
  invisible(x)
}

#' Find LD neighbors of a focal variant
#'
#' @param ld An `ld_result` object.
#' @param variant Variant ID, index, or genomic position.
#' @param threshold Minimum LD value.
#' @param metric `r2` or `dprime`.
#' @return Variant table ordered by decreasing LD.
#' @export
ld_neighbors <- function(ld, variant, threshold = 0.8, metric = c("r2", "dprime")) {
  if (!inherits(ld, "ld_result")) .stopf("ld must be an ld_result object.")
  metric <- match.arg(metric)
  mat <- ld[[metric]]
  if (is.null(mat)) .stopf("%s was not calculated.", metric)
  v <- ld$data$variants
  if (is.character(variant)) idx <- match(variant[1L], v$id)
  else if (length(variant) == 1L && variant %in% seq_len(nrow(v))) idx <- as.integer(variant)
  else idx <- which.min(abs(v$pos - as.numeric(variant[1L])))
  if (is.na(idx)) .stopf("Focal variant was not found.")
  keep <- is.finite(mat[idx, ]) & mat[idx, ] >= threshold
  out <- v[keep, , drop = FALSE]
  out[[metric]] <- mat[idx, keep]
  out$distance <- out$pos - v$pos[idx]
  out <- out[order(-out[[metric]], abs(out$distance)), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "focal_variant") <- v$id[idx]
  out
}

#' Calculate LD separately by population group
#'
#' @param x An `ld_data` object.
#' @param groups Named group vector, group vector in sample order, or a two-column sample/group data frame.
#' @param ... Arguments passed to `ld_compute`.
#' @return A named list of `ld_result` objects with class `ld_grouped`.
#' @export
ld_by_group <- function(x, groups, ...) {
  if (!inherits(x, "ld_data")) .stopf("x must be an ld_data object.")
  if (is.data.frame(groups)) {
    if (ncol(groups) < 2L) .stopf("groups data frame needs sample and group columns.")
    g <- setNames(as.character(groups[[2L]]), as.character(groups[[1L]]))
    groups <- g[x$samples]
  } else {
    groups <- as.character(groups)
    if (!is.null(names(groups))) groups <- groups[x$samples]
  }
  if (length(groups) != x$n_samples || anyNA(groups)) .stopf("Every sample must resolve to exactly one group.")
  lev <- unique(groups)
  if (length(lev) < 2L) {
    counts <- table(groups, useNA = "ifany")
    summary <- paste(sprintf("%s=%d", names(counts), as.integer(counts)), collapse = ", ")
    .stopf("At least two groups are required after matching x$samples; found %d group (%s). Read x without a single-group samples filter and ensure the group-file sample IDs overlap x$samples.",
           length(lev), summary)
  }
  out <- setNames(lapply(lev, function(g) ld_compute(subset_ld_data(x, samples = groups == g), ...)), lev)
  class(out) <- c("ld_grouped", "list")
  out
}

#' Compare LD matrices between groups
#'
#' @param grouped Result from `ld_by_group`.
#' @param metric `r2` or `dprime`.
#' @param reference Name or index of the reference group.
#' @return A list of group-minus-reference difference matrices.
#' @export
ld_compare <- function(grouped, metric = c("r2", "dprime"), reference = 1L) {
  if (!inherits(grouped, "ld_grouped")) .stopf("grouped must be returned by ld_by_group().")
  metric <- match.arg(metric)
  if (is.character(reference)) reference <- match(reference[1L], names(grouped))
  reference <- as.integer(reference[1L])
  if (is.na(reference) || reference < 1L || reference > length(grouped)) .stopf("Invalid reference group.")
  ref <- grouped[[reference]][[metric]]
  if (is.null(ref)) .stopf("%s is not available in the reference result.", metric)
  idx <- setdiff(seq_along(grouped), reference)
  out <- setNames(lapply(idx, function(i) {
    if (is.null(grouped[[i]][[metric]])) .stopf("%s is not available for group %s.", metric, names(grouped)[i])
    grouped[[i]][[metric]] - ref
  }), paste0(names(grouped)[idx], "_minus_", names(grouped)[reference]))
  class(out) <- c("ld_comparison", "list")
  attr(out, "metric") <- metric
  out
}

#' Estimate haplotype frequencies
#'
#' @param x An `ld_data` object.
#' @param variants Variant IDs, a single `|`-delimited block string, or indices;
#'   all variants by default.
#' @param min_frequency Minimum reported frequency.
#' @param drop_missing Drop chromosomes/samples with missing alleles.
#' @return A haplotype frequency table.
#' @export
haplotype_frequencies <- function(x, variants = NULL, min_frequency = 0.01,
                                  drop_missing = TRUE) {
  if (!inherits(x, "ld_data")) .stopf("x must be an ld_data object.")
  if (is.null(variants)) idx <- seq_len(x$n_variants)
  else if (is.character(variants)) {
    variants <- as.character(variants)
    if (length(variants) == 1L && !is.na(variants) && grepl("|", variants, fixed = TRUE)) {
      variants <- strsplit(variants, "|", fixed = TRUE)[[1L]]
      variants <- trimws(variants)
    }
    variants <- variants[!is.na(variants) & nzchar(variants)]
    idx <- match(variants, x$variants$id)
  }
  else idx <- as.integer(variants)
  if (!length(idx)) {
    .stopf("No variants were supplied. Check that at least one LD block was detected (nrow(blocks) > 0) before calling haplotype_frequencies().")
  }
  if (anyNA(idx) || any(idx < 1L | idx > x$n_variants)) .stopf("Invalid variants.")
  if (length(idx) > 30L) .warnf("Inferring frequencies across %d variants may produce many rare haplotypes.", length(idx))
  phase_ok <- !is.null(x$haplotypes) &&
    (!"phase_complete" %in% names(x$variants) || all(x$variants$phase_complete[idx]))
  if (phase_ok) {
    h <- x$haplotypes[, idx, drop = FALSE]
  } else {
    g <- x$genotypes[, idx, drop = FALSE]
    g[g == 1] <- NA_real_
    h <- g / 2
    .warnf("Complete phase was unavailable; frequencies use homozygous sample haplotypes only.")
  }
  if (drop_missing) h <- h[rowSums(is.na(h)) == 0L, , drop = FALSE]
  if (!nrow(h)) .stopf("No complete haplotypes remain.")
  strings <- apply(h, 1L, function(z) paste(ifelse(is.na(z), ".", z), collapse = ""))
  counts <- sort(table(strings), decreasing = TRUE)
  out <- data.frame(haplotype = names(counts), count = as.integer(counts),
                    frequency = as.numeric(counts) / sum(counts), stringsAsFactors = FALSE)
  out <- out[out$frequency >= min_frequency, , drop = FALSE]
  attr(out, "variants") <- x$variants$id[idx]
  attr(out, "n_chromosomes") <- length(strings)
  rownames(out) <- NULL
  out
}

#' Summarize LD decay with physical distance
#'
#' @param ld An `ld_result` object.
#' @param metric `r2` or `dprime`.
#' @param n_bins Number of equal-width distance bins.
#' @param bin_width Optional fixed bin width in base pairs.
#' @param max_distance Optional maximum distance.
#' @return A data frame of binned LD summaries.
#' @export
ld_decay <- function(ld, metric = c("r2", "dprime"), n_bins = 30L,
                     bin_width = NULL, max_distance = NULL) {
  if (!inherits(ld, "ld_result")) .stopf("ld must be an ld_result object.")
  metric <- match.arg(metric)
  mat <- ld[[metric]]
  if (is.null(mat)) .stopf("%s was not calculated.", metric)
  v <- ld$data$variants
  ij <- which(upper.tri(mat) & is.finite(mat), arr.ind = TRUE)
  if (!nrow(ij)) .stopf("No finite pairwise LD values are available.")
  distance <- abs(v$pos[ij[, 2L]] - v$pos[ij[, 1L]])
  value <- mat[ij]
  if (is.null(max_distance)) max_distance <- max(distance)
  keep <- distance <= max_distance
  distance <- distance[keep]
  value <- value[keep]
  if (is.null(bin_width)) bin_width <- max(1, max_distance / n_bins)
  bin <- floor(distance / bin_width)
  lev <- sort(unique(bin))
  out <- do.call(rbind, lapply(lev, function(b) {
    z <- value[bin == b]
    data.frame(distance_start = b * bin_width,
               distance_end = min((b + 1) * bin_width, max_distance),
               distance_mid = (b + 0.5) * bin_width,
               mean = mean(z, na.rm = TRUE), median = stats::median(z, na.rm = TRUE),
               q25 = stats::quantile(z, 0.25, na.rm = TRUE, names = FALSE),
               q75 = stats::quantile(z, 0.75, na.rm = TRUE, names = FALSE),
               n_pairs = length(z), stringsAsFactors = FALSE)
  }))
  attr(out, "metric") <- metric
  class(out) <- c("ld_decay", class(out))
  out
}

#' Plot an LD decay summary
#'
#' @param x Result from `ld_decay`.
#' @param col Line color.
#' @param ribbon_col Interquartile ribbon color.
#' @param xlab,ylab Axis labels.
#' @param ... Additional graphical parameters.
#' @return The input invisibly.
#' @export
plot_ld_decay <- function(x, col = "#B2182B", ribbon_col = "#FDD0A2",
                          xlab = "Physical distance (bp)", ylab = NULL, ...) {
  if (!inherits(x, "ld_decay")) .stopf("x must be returned by ld_decay().")
  metric <- attr(x, "metric")
  if (is.null(ylab)) ylab <- if (metric == "r2") expression(r^2) else expression(D * "'")
  graphics::plot(x$distance_mid, x$median, type = "n", xlab = xlab, ylab = ylab,
                 ylim = range(c(0, x$q25, x$q75), finite = TRUE), ...)
  graphics::polygon(c(x$distance_mid, rev(x$distance_mid)), c(x$q25, rev(x$q75)),
                    col = ribbon_col, border = NA)
  graphics::lines(x$distance_mid, x$median, col = col, lwd = 2)
  graphics::points(x$distance_mid, x$median, col = col, pch = 16, cex = 0.6)
  invisible(x)
}
