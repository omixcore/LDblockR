.prefix2 <- function(x) {
  n <- nrow(x)
  out <- matrix(0, n + 1L, n + 1L)
  for (i in seq_len(n)) out[i + 1L, -1L] <- out[i, -1L] + cumsum(x[i, ])
  out
}

.interval_pair_sum <- function(prefix, i, j) {
  (prefix[j + 1L, j + 1L] - prefix[i, j + 1L] -
     prefix[j + 1L, i] + prefix[i, i]) / 2
}

.candidate_table <- function(candidates, variants, method) {
  if (!length(candidates)) return(.empty_blocks())
  out <- do.call(rbind, lapply(candidates, function(z) {
    idx <- z$i:z$j
    data.frame(
      block = NA_character_, chr = as.character(variants$chr[z$i]),
      start = variants$pos[z$i], end = variants$pos[z$j],
      length_bp = variants$pos[z$j] - variants$pos[z$i] + 1,
      n_snps = length(idx), start_index = z$i, end_index = z$j,
      snps = paste(variants$id[idx], collapse = "|"), method = method,
      strong_fraction = .null_coalesce(z$fraction, NA_real_),
      stringsAsFactors = FALSE
    )
  }))
  out <- out[order(-out$length_bp, -out$n_snps, out$start), , drop = FALSE]
  occupied <- rep(FALSE, nrow(variants))
  keep <- logical(nrow(out))
  for (k in seq_len(nrow(out))) {
    idx <- out$start_index[k]:out$end_index[k]
    if (!any(occupied[idx])) {
      keep[k] <- TRUE
      occupied[idx] <- TRUE
    }
  }
  out <- out[keep, , drop = FALSE]
  out <- out[order(.chromosome_rank(out$chr), out$chr, out$start, out$end), , drop = FALSE]
  out$block <- paste0("Block", seq_len(nrow(out)))
  rownames(out) <- NULL
  class(out) <- c("ld_blocks", class(out))
  out
}

.enumerate_fraction_blocks <- function(strong, informative, positions,
                                       fraction_cut, max_span, min_snps,
                                       strict = FALSE) {
  n <- nrow(strong)
  diag(strong) <- FALSE
  diag(informative) <- FALSE
  ps <- .prefix2(strong * 1)
  pi <- .prefix2(informative * 1)
  candidates <- list()
  k <- 0L
  for (i in seq_len(n - min_snps + 1L)) {
    for (j in seq.int(i + min_snps - 1L, n)) {
      if (positions[j] - positions[i] > max_span) break
      ni <- .interval_pair_sum(pi, i, j)
      if (ni <= 0) next
      ns <- .interval_pair_sum(ps, i, j)
      ok <- if (strict) ns > fraction_cut * ni else ns / ni >= fraction_cut
      if (ok) {
        k <- k + 1L
        candidates[[k]] <- list(i = i, j = j, fraction = ns / ni)
      }
    }
  }
  candidates
}

.four_gamete_matrix <- function(x, indices, min_haplotype_frequency) {
  n <- length(indices)
  recomb <- matrix(FALSE, n, n)
  informative <- matrix(FALSE, n, n)
  pc <- if ("phase_complete" %in% names(x$variants)) x$variants$phase_complete else rep(FALSE, x$n_variants)
  for (a in seq_len(n - 1L)) {
    i <- indices[a]
    for (b in (a + 1L):n) {
      j <- indices[b]
      use_h <- !is.null(x$haplotypes) && isTRUE(pc[i]) && isTRUE(pc[j])
      s <- if (use_h) .pair_haplotype_stats(h1 = x$haplotypes[, i], h2 = x$haplotypes[, j]) else
        .pair_haplotype_stats(x$genotypes[, i], x$genotypes[, j])
      if (all(is.finite(s$freq))) {
        informative[a, b] <- informative[b, a] <- TRUE
        recomb[a, b] <- recomb[b, a] <- all(s$freq >= min_haplotype_frequency)
      }
    }
  }
  list(recombination = recomb, informative = informative)
}

#' Detect haplotype blocks
#'
#' @param ld An `ld_result` object.
#' @param method `gabriel`, `solid_spine`, `strong`, `four_gamete`, `fixed`, or `none`.
#' @param metric LD metric for non-Gabriel methods.
#' @param strong_cut Strong-LD threshold for the `strong` method.
#' @param strong_fraction Required strong-pair fraction.
#' @param spine_cut Endpoint LD threshold for solid spine.
#' @param max_span Maximum block length in base pairs.
#' @param min_snps Minimum SNPs per block.
#' @param fixed Data frame of chromosome, start, and end coordinates.
#' @param gabriel_lowci,gabriel_highci,gabriel_recomb_highci Gabriel pair-class thresholds.
#' @param gabriel_inform_fraction Required strong fraction among informative pairs.
#' @param ci_step Likelihood grid step when D-prime intervals must be calculated.
#' @param four_gamete_min_frequency Minimum frequency for each of four haplotypes.
#' @return A data frame of non-overlapping blocks with class `ld_blocks`.
#' @export
detect_ld_blocks <- function(ld,
                             method = c("gabriel", "solid_spine", "strong", "four_gamete", "fixed", "none"),
                             metric = c("dprime", "r2"),
                             strong_cut = 0.85, strong_fraction = 0.90,
                             spine_cut = 0.80, max_span = Inf, min_snps = 2L,
                             fixed = NULL,
                             gabriel_lowci = 0.70, gabriel_highci = 0.98,
                             gabriel_recomb_highci = 0.90,
                             gabriel_inform_fraction = 0.95,
                             ci_step = 0.01,
                             four_gamete_min_frequency = 0.01) {
  if (!inherits(ld, "ld_result")) .stopf("ld must be an ld_result object.")
  method <- match.arg(method)
  metric <- match.arg(metric)
  if (method == "none") return(.empty_blocks())
  v <- ld$data$variants
  if (method == "fixed") {
    if (is.null(fixed)) .stopf("fixed block coordinates are required for method='fixed'.")
    fixed <- as.data.frame(fixed, stringsAsFactors = FALSE)
    chr_col <- .match_column(fixed, c("chr", "chrom", "chromosome"), TRUE, "chromosome")
    start_col <- .match_column(fixed, c("start", "bp1", "from"), TRUE, "start")
    end_col <- .match_column(fixed, c("end", "bp2", "to"), TRUE, "end")
    cand <- list()
    k <- 0L
    for (r in seq_len(nrow(fixed))) {
      idx <- which(v$chr == fixed[[chr_col]][r] & v$pos >= fixed[[start_col]][r] & v$pos <= fixed[[end_col]][r])
      if (length(idx) >= min_snps) {
        k <- k + 1L
        cand[[k]] <- list(i = min(idx), j = max(idx), fraction = NA_real_)
      }
    }
    return(.candidate_table(cand, v, method))
  }
  if (method == "gabriel" && (is.null(ld$dprime) || is.null(ld$ci_lower))) {
    ld_ci <- ld_compute(ld$data, measure = "dprime", max_distance = ld$max_distance,
                        min_n = ld$min_n, ci = TRUE, ci_step = ci_step)
  } else {
    ld_ci <- ld
  }
  all_candidates <- list()
  candidate_count <- 0L
  chromosomes <- unique(as.character(v$chr))
  for (chr in chromosomes) {
    global_idx <- which(as.character(v$chr) == chr)
    if (length(global_idx) < min_snps) next
    pos <- v$pos[global_idx]
    local <- list()
    if (method == "gabriel") {
      low <- ld_ci$ci_lower[global_idx, global_idx, drop = FALSE]
      high <- ld_ci$ci_upper[global_idx, global_idx, drop = FALSE]
      strong <- is.finite(low) & is.finite(high) & low >= gabriel_lowci & high >= gabriel_highci
      recomb <- is.finite(high) & high < gabriel_recomb_highci
      informative <- strong | recomb
      local <- .enumerate_fraction_blocks(strong, informative, pos,
                                          gabriel_inform_fraction, max_span,
                                          min_snps, strict = TRUE)
    } else if (method == "strong") {
      mat <- ld[[metric]]
      if (is.null(mat)) .stopf("%s was not calculated.", metric)
      mat <- mat[global_idx, global_idx, drop = FALSE]
      strong <- is.finite(mat) & mat >= strong_cut
      informative <- is.finite(mat)
      local <- .enumerate_fraction_blocks(strong, informative, pos,
                                          strong_fraction, max_span, min_snps)
    } else if (method == "solid_spine") {
      mat <- ld[[metric]]
      if (is.null(mat)) .stopf("%s was not calculated.", metric)
      mat <- mat[global_idx, global_idx, drop = FALSE]
      k <- 0L
      for (i in seq_len(nrow(mat) - min_snps + 1L)) {
        for (j in seq.int(i + min_snps - 1L, nrow(mat))) {
          if (pos[j] - pos[i] > max_span) break
          left <- mat[i, i:j]
          right <- mat[i:j, j]
          if (all(is.finite(c(left, right))) && all(c(left, right) >= spine_cut)) {
            k <- k + 1L
            local[[k]] <- list(i = i, j = j,
                               fraction = mean(mat[i:j, i:j][upper.tri(mat[i:j, i:j])] >= spine_cut, na.rm = TRUE))
          }
        }
      }
    } else if (method == "four_gamete") {
      fg <- .four_gamete_matrix(ld$data, global_idx, four_gamete_min_frequency)
      strong <- fg$informative & !fg$recombination
      local <- .enumerate_fraction_blocks(strong, fg$informative, pos,
                                          1, max_span, min_snps)
    }
    if (length(local)) {
      for (z in local) {
        candidate_count <- candidate_count + 1L
        all_candidates[[candidate_count]] <- list(
          i = global_idx[z$i], j = global_idx[z$j], fraction = z$fraction
        )
      }
    }
  }
  .candidate_table(all_candidates, v, method)
}

#' @export
print.ld_blocks <- function(x, ...) {
  cat("<ld_blocks>", nrow(x), "non-overlapping block(s)\n")
  if (nrow(x)) {
    cat("  method:", paste(unique(x$method), collapse = ", "), "\n")
    cat("  variants in blocks:", sum(x$n_snps), "\n")
  }
  invisible(x)
}

#' Select tag SNPs by greedy LD coverage
#'
#' @param ld An `ld_result` object.
#' @param threshold Minimum LD for coverage.
#' @param metric `r2` or `dprime`.
#' @param blocks Optional block table; selection is then performed per block.
#'   Both `ld_blocks` output with `start_index`/`end_index` and coordinate
#'   tables with `chr`/`start`/`end` are accepted.
#' @param include_unblocked Select unblocked variants as singleton tags.
#' @return A tag-SNP table including covered variant IDs.
#' @export
select_tag_snps <- function(ld, threshold = 0.8, metric = c("r2", "dprime"),
                            blocks = NULL, include_unblocked = TRUE) {
  if (!inherits(ld, "ld_result")) .stopf("ld must be an ld_result object.")
  metric <- match.arg(metric)
  mat <- ld[[metric]]
  if (is.null(mat)) .stopf("%s was not calculated.", metric)
  v <- ld$data$variants
  if (is.null(blocks)) {
    sets <- list(All = seq_len(nrow(v)))
  } else {
    blocks <- as.data.frame(blocks, stringsAsFactors = FALSE)
    if (!nrow(blocks)) {
      sets <- list(All = seq_len(nrow(v)))
    } else {
      block_col <- .match_column(blocks, c("block", "name", "id"), FALSE)
      labels <- if (is.null(block_col)) paste0("Block", seq_len(nrow(blocks))) else
        as.character(blocks[[block_col]])
      labels[is.na(labels) | !nzchar(labels)] <- paste0("Block", which(is.na(labels) | !nzchar(labels)))
      if (all(c("start_index", "end_index") %in% names(blocks))) {
        sets <- lapply(seq_len(nrow(blocks)), function(i) {
          si <- suppressWarnings(as.integer(blocks$start_index[i]))
          ei <- suppressWarnings(as.integer(blocks$end_index[i]))
          if (any(!is.finite(c(si, ei))) || si < 1L || ei < si || ei > nrow(v)) {
            integer()
          } else {
            seq.int(si, ei)
          }
        })
      } else if (all(c("start", "end") %in% names(blocks))) {
        start_col <- .match_column(blocks, c("start", "bp1", "from"), TRUE, "block start")
        end_col <- .match_column(blocks, c("end", "bp2", "to"), TRUE, "block end")
        chr_col <- .match_column(blocks, c("chr", "#chr", "chrom", "chromosome"), FALSE)
        if (is.null(chr_col) && length(unique(as.character(v$chr))) != 1L) {
          .stopf("Coordinate block tables without a chromosome column are only valid for a single-chromosome LD region.")
        }
        v_chr <- as.character(v$chr)
        v_chr_key <- sub("^chr", "", v_chr, ignore.case = TRUE)
        sets <- lapply(seq_len(nrow(blocks)), function(i) {
          start <- suppressWarnings(as.numeric(blocks[[start_col]][i]))
          end <- suppressWarnings(as.numeric(blocks[[end_col]][i]))
          chr <- if (is.null(chr_col)) v_chr[1L] else as.character(blocks[[chr_col]][i])
          chr_key <- sub("^chr", "", chr, ignore.case = TRUE)
          if (any(!is.finite(c(start, end))) || end < start) return(integer())
          which(v_chr_key == chr_key & v$pos >= start & v$pos <= end)
        })
      } else {
        .stopf("blocks must contain start_index/end_index or chr/start/end columns.")
      }
      keep <- lengths(sets) > 0L
      if (!any(keep)) {
        .stopf("No LD variants fall inside the supplied block coordinates.")
      }
      sets <- sets[keep]
      labels <- labels[keep]
      names(sets) <- make.unique(labels)
    }
    covered_by_blocks <- unique(unlist(sets, use.names = FALSE))
    if (include_unblocked) {
      rest <- setdiff(seq_len(nrow(v)), covered_by_blocks)
      if (length(rest)) sets <- c(sets, setNames(as.list(rest), paste0("Unblocked_", v$id[rest])))
    }
  }
  result <- list()
  k <- 0L
  for (set_name in names(sets)) {
    idx <- sets[[set_name]]
    sub <- mat[idx, idx, drop = FALSE]
    cover <- is.finite(sub) & sub >= threshold
    diag(cover) <- TRUE
    uncovered <- rep(TRUE, length(idx))
    while (any(uncovered)) {
      scores <- colSums(cover[uncovered, , drop = FALSE])
      scores[!uncovered] <- -1
      best_score <- max(scores)
      choices <- which(scores == best_score)
      if (length(choices) > 1L) {
        maf <- v$maf[idx[choices]]
        call_count <- v$call_count[idx[choices]]
        choices <- choices[order(-maf, -call_count, v$pos[idx[choices]])]
      }
      tag_local <- choices[1L]
      newly <- which(uncovered & cover[, tag_local])
      uncovered[newly] <- FALSE
      tag_global <- idx[tag_local]
      k <- k + 1L
      result[[k]] <- data.frame(
        tag = v$id[tag_global], chr = v$chr[tag_global], pos = v$pos[tag_global],
        maf = v$maf[tag_global], block = set_name,
        n_covered = length(newly), covered = paste(v$id[idx[newly]], collapse = "|"),
        threshold = threshold, metric = metric, stringsAsFactors = FALSE
      )
    }
  }
  out <- if (length(result)) do.call(rbind, result) else data.frame()
  rownames(out) <- NULL
  class(out) <- c("ld_tags", class(out))
  out
}
