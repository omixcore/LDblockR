# Publication-ready association diagnostics

.gwas_plot_table <- function(gwas) {
  tab <- read_gwas(gwas)
  if (!nrow(tab)) .stopf("GWAS table has no finite association rows.")
  tab$chr <- as.character(tab$chr)
  tab$pos <- as.numeric(tab$pos)
  tab$p <- as.numeric(tab$p)
  tab$logp <- as.numeric(tab$logp)
  tab
}

.gwas_plot_highlight <- function(highlight, tab) {
  marked <- rep(FALSE, nrow(tab))
  if (is.null(highlight) || !length(highlight)) return(marked)
  if (is.data.frame(highlight)) {
    h <- as.data.frame(highlight, stringsAsFactors = FALSE)
    id_col <- .match_column(h, c("id", "snp", "marker", "rsid", "name"), FALSE)
    if (!is.null(id_col)) marked <- marked | as.character(tab$id) %in% as.character(h[[id_col]])
    chr_col <- .match_column(h, c("chr", "chrom", "chromosome"), FALSE)
    pos_col <- .match_column(h, c("pos", "position", "bp", "site"), FALSE)
    if (!is.null(chr_col) && !is.null(pos_col)) {
      key <- paste(tab$chr, tab$pos, sep = "\r")
      hkey <- paste(as.character(h[[chr_col]]), as.numeric(h[[pos_col]]), sep = "\r")
      marked <- marked | key %in% hkey
    }
    return(marked)
  }
  highlight <- as.character(highlight)
  marked <- marked | as.character(tab$id) %in% highlight
  if (any(grepl(":", highlight, fixed = TRUE))) {
    key <- paste(tab$chr, tab$pos, sep = ":")
    marked <- marked | key %in% highlight
  }
  marked
}

.gwas_point_sizes <- function(tab, point_size, point_size_by) {
  pve_col <- .match_column(tab, c("PVE", "pve", "var_explained",
                                  "variance_explained"), FALSE)
  pve <- if (is.null(pve_col)) rep(NA_real_, nrow(tab)) else
    suppressWarnings(as.numeric(tab[[pve_col]]))
  if (point_size_by != "PVE" || !any(is.finite(pve) & pve > 0)) {
    return(list(size = rep(point_size, nrow(tab)), pve = pve, max = NA_real_))
  }
  max_pve <- max(pve, na.rm = TRUE)
  scaled <- pmin(1, pmax(0, pve / max_pve))
  size <- point_size * (0.78 + 1.22 * sqrt(scaled))
  size[!is.finite(size)] <- point_size
  list(size = size, pve = pve, max = max_pve)
}

.gwas_thresholds <- function(cutline, threshold, n) {
  if (!is.null(threshold)) {
    if (!is.null(cutline) && length(cutline) && !all(is.na(cutline))) {
      .warnf("threshold overrides cutline.")
    }
    cutline <- threshold
  }
  if (is.null(cutline) || !length(cutline) || all(is.na(cutline))) {
    return(-log10(0.05 / n))
  }
  z <- suppressWarnings(as.numeric(cutline))
  z <- z[!is.na(z)]
  if (!length(z)) return(-log10(0.05 / n))
  if (any(!is.finite(z) & !is.infinite(z))) {
    .stopf("threshold values must be finite, Inf, or NULL.")
  }
  if (any(z < 0)) .stopf("threshold values must be non-negative.")
  z
}

.manhattan_chr_colors <- function(chr, chr_colors, palette, color_cycle = TRUE) {
  cols <- if (!is.null(chr_colors)) chr_colors else palette
  cols <- as.character(cols)
  if (!length(cols) || anyNA(cols) || any(!nzchar(cols))) {
    .stopf("chr_colors/palette must contain at least one valid color.")
  }
  names0 <- names(cols)
  has_names <- !is.null(names0) && any(!is.na(names0) & nzchar(names0))
  if (has_names) {
    named <- !is.na(names0) & nzchar(names0)
    by_chr <- cols[match(chr, names0)]
    missing <- is.na(by_chr)
    if (any(missing)) {
      fallback <- unname(cols[!named])
      if (!color_cycle && !length(fallback)) {
        .stopf("chr_colors does not define all chromosomes and color_cycle=FALSE.")
      }
      if (!color_cycle && length(fallback) < sum(missing)) {
        .stopf("chr_colors does not define all chromosomes and color_cycle=FALSE.")
      }
      if (!length(fallback)) fallback <- unname(cols)
      by_chr[missing] <- rep(fallback, length.out = sum(missing))
    }
  } else {
    if (!color_cycle && length(cols) < length(chr)) {
      .stopf("chr_colors/palette has %d colors for %d chromosomes; set color_cycle=TRUE or provide one color per chromosome.",
             length(cols), length(chr))
    }
    by_chr <- if (color_cycle) rep(unname(cols), length.out = length(chr)) else unname(cols)[seq_along(chr)]
  }
  stats::setNames(by_chr, chr)
}

.manhattan_chr_labels <- function(chr, labels) {
  if (is.null(labels)) return(chr)
  labels <- as.character(labels)
  nms <- names(labels)
  if (!is.null(nms) && any(!is.na(nms) & nzchar(nms))) {
    out <- labels[match(chr, nms)]
    out[is.na(out)] <- chr[is.na(out)]
    return(unname(out))
  }
  if (length(labels) != length(chr)) {
    .stopf("chr_labels must have one label per chromosome (%d) when unnamed; supply a named vector keyed by the chromosomes or subset gwas before plotting.",
           length(chr))
  }
  labels
}

.manhattan_colors <- function(tab, color_by, thresholds, palette,
                              chr_colors = NULL, color_cycle = TRUE,
                              color_significant = FALSE,
                              significant_color = "#E64B35",
                              nonsignificant_color = "#2C7FB8",
                              alpha = 0.82) {
  chr <- unique(tab$chr)
  if (color_by == "chromosome") {
    by_chr <- .manhattan_chr_colors(chr, chr_colors, palette, color_cycle)
    out <- unname(by_chr[match(tab$chr, chr)])
    if (isTRUE(color_significant) && is.finite(thresholds[1L])) {
      out[tab$logp >= thresholds[1L]] <- significant_color
    }
    return(grDevices::adjustcolor(out, alpha.f = alpha))
  }
  if (color_by == "significance") {
    out <- if (is.finite(thresholds[1L])) {
      ifelse(tab$logp >= thresholds[1L], significant_color, nonsignificant_color)
    } else rep(nonsignificant_color, nrow(tab))
    return(grDevices::adjustcolor(out, alpha.f = alpha))
  }
  pve_col <- .match_column(tab, c("PVE", "pve", "var_explained",
                                  "variance_explained"), FALSE)
  pve <- if (is.null(pve_col)) rep(NA_real_, nrow(tab)) else
    suppressWarnings(as.numeric(tab[[pve_col]]))
  if (!any(is.finite(pve))) return(grDevices::adjustcolor(rep(nonsignificant_color, nrow(tab)), alpha.f = alpha))
  top <- max(pve, na.rm = TRUE)
  z <- pmin(1, pmax(0, pve / max(top, 1e-12)))
  out <- .color_values(z, .ld_palette(c("#2C7FB8", "#91CF60", "#D73027")))
  out[!is.finite(pve)] <- "#BDBDBD"
  grDevices::adjustcolor(out, alpha.f = alpha)
}

.manhattan_chr_lengths <- function(values, chr, observed, name) {
  if (is.null(values)) return(NULL)
  values <- suppressWarnings(as.numeric(values))
  if (anyNA(values) || any(!is.finite(values)) || any(values <= 0)) {
    .stopf("%s must contain positive finite chromosome lengths.", name)
  }
  nms <- names(values)
  has_names <- !is.null(nms) && any(!is.na(nms) & nzchar(nms))
  if (has_names) {
    out <- values[match(chr, nms)]
    if (anyNA(out)) {
      .stopf("%s must define every chromosome when named.", name)
    }
  } else {
    if (length(values) != length(chr)) {
      .stopf("%s must have one value per chromosome (%d) when unnamed; supply a named vector for all chromosomes or subset gwas before plotting.",
             name, length(chr))
    }
    out <- values
  }
  if (any(out < observed)) {
    bad <- paste(chr[out < observed], collapse = ", ")
    .warnf("%s is shorter than the largest observed coordinate for chromosome(s): %s; using the observed coordinate as the minimum width.",
           name, bad)
    out <- pmax(out, observed)
  }
  unname(out)
}

.prepare_manhattan <- function(tab, gap_ratio = 0.05, gap_bp = NULL,
                               gap_mode = c("fixed", "relative"),
                               chr_lengths = NULL,
                               width_mode = c("length", "observed", "equal")) {
  gap_mode <- match.arg(gap_mode)
  width_mode <- match.arg(width_mode)
  tab <- tab[is.finite(tab$pos) & is.finite(tab$logp), , drop = FALSE]
  if (!nrow(tab)) .stopf("No finite GWAS positions and P values remain.")
  tab <- tab[order(.chromosome_rank(tab$chr), tab$chr, tab$pos, tab$id), , drop = FALSE]
  rownames(tab) <- NULL
  chr <- unique(tab$chr)
  ranges <- lapply(chr, function(z) range(tab$pos[tab$chr == z], finite = TRUE))
  observed_spans <- vapply(ranges, diff, numeric(1L))
  observed_spans[!is.finite(observed_spans) | observed_spans <= 0] <- 1
  observed_extent <- vapply(ranges, function(z) max(z[2L], 1), numeric(1L))
  spans <- .manhattan_chr_lengths(chr_lengths, chr, observed_extent, "chr_lengths")
  if (is.null(spans)) {
    if (width_mode == "length") {
      # In a whole-genome table, the largest observed coordinate is the best
      # available proxy for the chromosome's physical length. The optional
      # chr_lengths argument can replace this with assembly lengths.
      spans <- observed_extent
      spans <- pmax(spans, observed_spans)
    } else if (width_mode == "equal") {
      spans <- rep(max(observed_spans), length(chr))
    } else {
      spans <- observed_spans
    }
  }
  gap_ratio <- as.numeric(gap_ratio[1L])
  if (!is.finite(gap_ratio) || gap_ratio < 0) .stopf("gap_ratio must be finite and non-negative.")
  if (!is.null(gap_bp)) {
    gap <- as.numeric(gap_bp[1L])
    if (!is.finite(gap) || gap < 0) .stopf("gap_bp must be finite and non-negative.")
    if (gap_mode != "fixed") {
      .warnf("gap_bp is supplied; using gap_mode='fixed'.")
      gap_mode <- "fixed"
    }
  } else if (gap_mode == "fixed") {
    # A single gap value is used at every chromosome boundary.  This is the
    # publication default and produces equal visual spacing even when the
    # chromosome widths are unequal.
    gap <- max(1, stats::median(spans) * gap_ratio)
  } else {
    gap <- pmax(1, spans[-length(spans)] * gap_ratio)
  }
  n_chr <- length(spans)
  gap_values <- if (n_chr <= 1L) numeric() else if (length(gap) == 1L) {
    rep(gap, n_chr - 1L)
  } else {
    gap
  }
  offsets <- if (n_chr <= 1L) 0 else
    c(0, cumsum(spans[-n_chr] + gap_values))
  names(offsets) <- chr
  # CMplot-style cumulative coordinates: when physical chromosome lengths
  # determine the panel, retain the original within-chromosome bp coordinate
  # and add the cumulative offset. Only observed-span/equal layouts are
  # re-based to the first observed marker.
  origins <- vapply(ranges, function(z) z[1L], numeric(1L))
  if (width_mode == "length" || !is.null(chr_lengths)) origins[] <- 0
  names(origins) <- chr
  tab$x <- tab$pos - origins[tab$chr] + offsets[tab$chr]
  mids <- offsets + spans / 2
  names(mids) <- chr
  list(data = tab, chr = chr, ranges = ranges, spans = spans,
       observed_spans = observed_spans, offsets = offsets, mids = mids,
       gap = gap, gap_values = gap_values, gap_mode = gap_mode,
       width_mode = width_mode, chr_lengths = spans)
}

.draw_manhattan <- function(x) {
  tab <- x$data
  threshold <- as.numeric(x$threshold)
  suggestive <- as.numeric(x$suggestive)
  finite_thresholds <- c(threshold[is.finite(threshold)],
                         suggestive[is.finite(suggestive)])
  y_max <- max(c(tab$logp, finite_thresholds), na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  y_lim <- if (is.null(x$ylim)) c(0, max(1, y_max * 1.10)) else x$ylim
  if (is.null(x$xlim)) {
    left <- x$chr_mid - x$chr_spans / 2
    right <- x$chr_mid + x$chr_spans / 2
    panel_lim <- range(c(left, right), finite = TRUE)
    pad <- max(diff(panel_lim) * 0.012, 1)
    x_lim <- panel_lim + c(-pad, pad)
  } else x_lim <- x$xlim
  if (diff(x_lim) <= 0) x_lim <- x_lim + c(-0.5, 0.5)
  graphics::plot(NA, NA, xlim = x_lim, ylim = y_lim, xaxs = "i", yaxs = "i",
                 axes = FALSE, xlab = "", ylab = "")
  if (isTRUE(x$show_grid)) {
    grid_at <- pretty(y_lim, n = 5L)
    grid_at <- grid_at[grid_at > 0 & grid_at < y_lim[2L]]
    if (length(grid_at)) graphics::abline(h = grid_at, col = x$grid_color,
                                           lty = x$grid_lty, lwd = x$grid_lwd)
  }
  graphics::abline(h = 0, col = "#666666", lwd = 0.8)
  signal <- if (length(threshold) && is.finite(threshold[1L])) {
    tab$logp >= threshold[1L]
  } else rep(FALSE, nrow(tab))
  point_cex <- x$point_size * ifelse(signal, x$signal_cex, 1)
  graphics::points(tab$x, tab$logp, pch = x$pch, cex = point_cex,
                   col = x$point_colors)
  if (any(x$highlight)) {
    graphics::points(tab$x[x$highlight], tab$logp[x$highlight], pch = x$highlight_pch,
                     cex = point_cex[x$highlight] * x$highlight_cex,
                     bg = x$highlight_col, col = "white", lwd = 0.7)
  }
  if (length(suggestive)) {
    for (i in seq_along(suggestive)) {
      if (is.finite(suggestive[i])) graphics::abline(
        h = suggestive[i],
        lty = rep(x$suggestive_lty, length.out = length(suggestive))[i],
        col = rep(x$suggestive_col, length.out = length(suggestive))[i],
        lwd = rep(x$suggestive_lwd, length.out = length(suggestive))[i])
    }
  }
  if (length(threshold)) {
    for (i in seq_along(threshold)) {
      if (is.finite(threshold[i])) graphics::abline(
        h = threshold[i],
        lty = rep(x$threshold_lty, length.out = length(threshold))[i],
        col = rep(x$threshold_col, length.out = length(threshold))[i],
        lwd = rep(x$threshold_lwd, length.out = length(threshold))[i])
    }
  }
  graphics::axis(1, at = x$chr_mid, labels = x$chr_labels, las = x$x_axis_las,
                 cex.axis = x$axis_cex,
                 col.axis = "#333333", col.ticks = "#555555", tck = -0.018)
  graphics::axis(2, las = 1, cex.axis = x$axis_cex)
  graphics::mtext(x$xlab, side = 1, line = 1.55, cex = x$label_cex)
  graphics::mtext(x$ylab, side = 2, line = 2.45, cex = x$label_cex)
  if (x$label_top > 0L) {
    ord <- order(tab$logp, decreasing = TRUE)
    ord <- ord[seq_len(min(x$label_top, length(ord)))]
    graphics::text(tab$x[ord], tab$logp[ord] + x$label_offset * diff(y_lim),
                   labels = tab$id[ord], cex = x$label_cex, srt = x$label_angle,
                   adj = c(0, 0), xpd = NA)
  }
  if (isTRUE(x$show_pve_legend) && is.finite(x$pve_max) && x$pve_max > 0) {
    lv <- c(0, x$pve_max / 2, x$pve_max)
    graphics::legend("topright", legend = formatC(lv, digits = 2, format = "f"),
                     title = "PVE", pch = 16,
                     pt.cex = x$base_point_size * (0.78 + 1.22 * sqrt(lv / x$pve_max)),
                     col = "#4D4D4D", bty = "n", cex = x$legend_cex,
                     inset = c(0.01, 0.01))
  }
  if (isTRUE(x$show_threshold_legend)) {
    threshold_ok <- is.finite(threshold)
    suggestive_ok <- is.finite(suggestive)
    n_lines <- sum(threshold_ok) + sum(suggestive_ok)
    if (n_lines) {
      labels <- c(x$threshold_labels[threshold_ok],
                  x$suggestive_labels[suggestive_ok])
      lty <- c(rep(x$threshold_lty, length.out = length(threshold))[threshold_ok],
               rep(x$suggestive_lty, length.out = length(suggestive))[suggestive_ok])
      col <- c(rep(x$threshold_col, length.out = length(threshold))[threshold_ok],
               rep(x$suggestive_col, length.out = length(suggestive))[suggestive_ok])
      lwd <- c(rep(x$threshold_lwd, length.out = length(threshold))[threshold_ok],
               rep(x$suggestive_lwd, length.out = length(suggestive))[suggestive_ok])
      graphics::legend("topleft", legend = labels, lty = lty, col = col,
                       lwd = lwd, bty = "n", cex = x$legend_cex,
                       inset = c(0.01, 0.01))
    }
  }
  if (!is.null(x$title) && nzchar(x$title)) {
    graphics::mtext(x$title, side = 3, line = 0.65, cex = x$title_cex, font = 2)
  }
  graphics::box(col = "#666666")
  invisible(x)
}

#' Draw a publication-ready Manhattan plot
#'
#' @param gwas A GWAS table, a path accepted by [read_gwas()], or the result of
#'   [gwas_mlm()]. Required columns are chromosome, position, and P/value.
#' @param cutline A `-log10(P)` threshold. `NULL` uses the Bonferroni threshold
#'   `-log10(0.05/n)`, and `Inf` suppresses the line. Kept for backward
#'   compatibility; use `threshold` for one or more lines.
#' @param suggestive Optional vector of suggestive `-log10(P)` thresholds.
#' @param palette Fallback chromosome colors. Values are recycled when
#'   `color_cycle = TRUE`.
#' @param color_by Point colors: chromosome, significance, or PVE.
#' @param point_size_by Encode PVE in point size when available, or use none.
#' @param point_size Base point size.
#' @param highlight SNP IDs, `chr:position` strings, or a table with ID/position.
#' @param highlight_col Color used for highlighted SNPs.
#' @param label_top Number of top SNPs to label.
#' @param show_pve_legend Add a PVE size legend when PVE is available.
#' @param title,xlab,ylab Plot title and axis labels.
#' @param draw Draw immediately; otherwise return an object for later printing.
#' @param threshold One or more `-log10(P)` thresholds. Overrides `cutline`.
#' @param threshold_col,threshold_color Colors for threshold lines; the latter
#'   is a compatibility alias. Values are recycled across thresholds.
#' @param threshold_lty,threshold_lwd Line types and widths for threshold lines.
#' @param suggestive_col,suggestive_color Colors for suggestive lines.
#' @param suggestive_lty,suggestive_lwd Line types and widths for suggestive lines.
#' @param chr_colors A named or unnamed vector of per-chromosome colors.
#'   Named colors are matched to chromosome labels.
#' @param color_cycle Recycle `chr_colors` or `palette` over chromosomes.
#' @param chr_labels Optional named or positional chromosome labels.
#' @param colors Alias for `chr_colors`.
#' @param gap_ratio Scaling factor used to derive an automatic gap when
#'   `gap_bp` is not supplied. With the default `gap_mode = "fixed"`, the
#'   resulting single gap is used at every chromosome boundary.
#' @param gap_bp Explicit fixed gap on the genomic-position scale between
#'   chromosomes. This is in the same coordinate unit as `pos` (for example,
#'   5e6 for a 5-Mb gap).
#' @param gap_mode Gap rule: `"fixed"` uses one identical gap at every
#'   chromosome boundary (the publication default); `"relative"` scales each
#'   boundary gap by the width of the preceding chromosome. Supplying `gap_bp`
#'   always selects the fixed rule.
#' @param chr_lengths Optional named or positional physical chromosome lengths
#'   in bp. They determine the relative panel widths and must be at least as
#'   large as the observed marker span.
#' @param width_mode Width rule when `chr_lengths` is `NULL`: `"length"`
#'   uses the largest observed coordinate, `"observed"` uses the marker span,
#'   and `"equal"` gives every chromosome the same width.
#' @param ylim,xlim Optional two-element plotting ranges.
#' @param pch Point plotting symbol.
#' @param point_cex Alias for `point_size`.
#' @param signal_cex Multiplicative size for points at or above the first threshold.
#' @param color_significant Recolor significant chromosome points.
#' @param significant_color,nonsignificant_color Colors used by significance mode.
#' @param alpha,alpha_points Point transparency in `[0, 1]`; the latter is an alias.
#' @param line_width Set both threshold and suggestive line widths.
#' @param threshold_labels,suggestive_labels Legend labels for the lines.
#' @param show_threshold_legend Add a line legend.
#' @param show_grid,grid_color,grid_lty,grid_lwd Horizontal grid-line controls.
#' @param axis_cex,label_cex,title_cex Font scaling controls.
#' @param label_angle,label_offset Controls for top-SNP labels.
#' @param x_axis_las Orientation of chromosome-axis labels.
#' @param highlight_pch,highlight_cex Symbol and size multiplier for highlights.
#' @param legend_cex Legend font scaling.
#' @return An object of class `manhattan_plot`.
#' @export
plot_manhattan <- function(gwas, cutline = NULL, suggestive = NULL,
                           palette = c("#2C7FB8", "#7B3294"),
                           color_by = c("chromosome", "significance", "pve"),
                           point_size_by = c("PVE", "none"), point_size = 0.58,
                           highlight = NULL, highlight_col = "#D73027",
                           label_top = 0L, show_pve_legend = TRUE,
                           title = "Manhattan plot", xlab = "Chromosome",
                           ylab = expression(-log[10](P)), draw = TRUE,
                           threshold = NULL, threshold_col = NULL,
                           threshold_lty = 2, threshold_lwd = 1,
                           threshold_color = NULL, suggestive_col = "#8C8C8C",
                           suggestive_lty = 3, suggestive_lwd = 0.9,
                           suggestive_color = NULL, chr_colors = NULL,
                           color_cycle = TRUE, chr_labels = NULL, colors = NULL,
                           gap_ratio = 0.05, gap_bp = NULL,
                           gap_mode = c("fixed", "relative"), chr_lengths = NULL,
                           width_mode = c("length", "observed", "equal"),
                           ylim = NULL,
                           xlim = NULL, pch = 16, point_cex = NULL,
                           signal_cex = 1, color_significant = FALSE,
                           significant_color = "#E64B35",
                           nonsignificant_color = "#2C7FB8", alpha = 0.82,
                           alpha_points = NULL, line_width = NULL,
                           threshold_labels = NULL, suggestive_labels = NULL,
                           show_threshold_legend = TRUE, show_grid = TRUE,
                           grid_color = "#EEF2F5", grid_lty = 1, grid_lwd = 0.7,
                           axis_cex = 0.72, label_cex = 0.58, title_cex = 1.04,
                           label_angle = 35, label_offset = 0.025,
                           x_axis_las = 1, highlight_pch = 21,
                           highlight_cex = 1.25, legend_cex = 0.62) {
  color_by <- match.arg(color_by)
  point_size_by <- match.arg(point_size_by)
  label_top <- max(0L, as.integer(label_top))
  if (!is.null(colors)) {
    if (!is.null(chr_colors)) .stopf("Use either colors or chr_colors, not both.")
    chr_colors <- colors
  }
  if (!is.null(threshold_color)) threshold_col <- threshold_color
  if (is.null(threshold_col)) threshold_col <- "#D73027"
  if (!is.null(suggestive_color)) suggestive_col <- suggestive_color
  if (!is.null(line_width)) {
    threshold_lwd <- line_width
    suggestive_lwd <- line_width
  }
  if (!is.null(alpha_points)) alpha <- alpha_points
  if (!is.null(point_cex)) point_size <- point_cex
  scalar_positive <- function(z, nm) {
    z <- as.numeric(z[1L])
    if (!is.finite(z) || z <= 0) .stopf("%s must be a positive finite number.", nm)
    z
  }
  point_size <- scalar_positive(point_size, "point_size")
  signal_cex <- scalar_positive(signal_cex, "signal_cex")
  highlight_cex <- scalar_positive(highlight_cex, "highlight_cex")
  pch <- as.numeric(pch[1L])
  highlight_pch <- as.numeric(highlight_pch[1L])
  if (!is.finite(pch) || !is.finite(highlight_pch)) .stopf("pch values must be finite.")
  alpha <- as.numeric(alpha[1L])
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) .stopf("alpha must be between 0 and 1.")
  if (!is.logical(color_cycle) || length(color_cycle) != 1L || is.na(color_cycle)) {
    .stopf("color_cycle must be TRUE or FALSE.")
  }
  validate_range <- function(z, nm, allow_null = TRUE) {
    if (is.null(z) && allow_null) return(NULL)
    z <- as.numeric(z)
    if (length(z) != 2L || any(!is.finite(z)) || z[1L] >= z[2L]) {
      .stopf("%s must be a finite increasing two-element range.", nm)
    }
    z
  }
  ylim <- validate_range(ylim, "ylim")
  xlim <- validate_range(xlim, "xlim")
  gap_ratio <- as.numeric(gap_ratio[1L])
  if (!is.finite(gap_ratio) || gap_ratio < 0) .stopf("gap_ratio must be finite and non-negative.")
  if (!is.null(gap_bp)) {
    gap_bp <- as.numeric(gap_bp[1L])
    if (!is.finite(gap_bp) || gap_bp < 0) .stopf("gap_bp must be finite and non-negative.")
  }
  width_mode <- match.arg(width_mode)
  gap_mode <- match.arg(gap_mode)
  tab <- .gwas_plot_table(gwas)
  prep <- .prepare_manhattan(tab, gap_ratio = gap_ratio, gap_bp = gap_bp,
                             gap_mode = gap_mode, chr_lengths = chr_lengths,
                             width_mode = width_mode)
  tab <- prep$data
  threshold <- .gwas_thresholds(cutline, threshold, nrow(tab))
  if (is.null(suggestive) || !length(suggestive)) {
    suggestive <- numeric()
  } else {
    suggestive <- suppressWarnings(as.numeric(suggestive))
    if (anyNA(suggestive) || any(suggestive < 0) || any(is.nan(suggestive))) {
      .stopf("suggestive must contain non-negative numbers.")
    }
  }
  threshold_lty <- as.numeric(threshold_lty)
  threshold_lwd <- as.numeric(threshold_lwd)
  suggestive_lty <- as.numeric(suggestive_lty)
  suggestive_lwd <- as.numeric(suggestive_lwd)
  if (!length(threshold_lty) || any(!is.finite(threshold_lty)) ||
      !length(threshold_lwd) || any(!is.finite(threshold_lwd)) || any(threshold_lwd <= 0) ||
      !length(suggestive_lty) || any(!is.finite(suggestive_lty)) ||
      !length(suggestive_lwd) || any(!is.finite(suggestive_lwd)) || any(suggestive_lwd <= 0)) {
    .stopf("Line types must be finite and line widths must be positive.")
  }
  threshold_col <- as.character(threshold_col)
  suggestive_col <- as.character(suggestive_col)
  if (!length(threshold_col) || anyNA(threshold_col) || !length(suggestive_col) || anyNA(suggestive_col)) {
    .stopf("Line colors must be non-missing character values.")
  }
  if (is.null(threshold_labels)) {
    threshold_labels <- ifelse(seq_along(threshold) == 1L, "Bonferroni", paste0("Threshold ", seq_along(threshold)))
  } else threshold_labels <- as.character(threshold_labels)
  if (length(threshold_labels) != length(threshold)) {
    .stopf("threshold_labels must have one label per threshold.")
  }
  if (is.null(suggestive_labels)) {
    suggestive_labels <- if (length(suggestive)) paste0("Suggestive ", seq_along(suggestive)) else character()
  } else suggestive_labels <- as.character(suggestive_labels)
  if (length(suggestive_labels) != length(suggestive)) {
    .stopf("suggestive_labels must have one label per suggestive threshold.")
  }
  numeric_positive <- function(z, nm) {
    z <- as.numeric(z[1L])
    if (!is.finite(z) || z <= 0) .stopf("%s must be positive and finite.", nm)
    z
  }
  axis_cex <- numeric_positive(axis_cex, "axis_cex")
  label_cex <- numeric_positive(label_cex, "label_cex")
  title_cex <- numeric_positive(title_cex, "title_cex")
  legend_cex <- numeric_positive(legend_cex, "legend_cex")
  label_offset <- as.numeric(label_offset[1L])
  if (!is.finite(label_offset) || label_offset < 0) .stopf("label_offset must be non-negative.")
  prep_chr_labels <- .manhattan_chr_labels(prep$chr, chr_labels)
  sizes <- .gwas_point_sizes(tab, point_size, point_size_by)
  point_colors <- .manhattan_colors(tab, color_by, threshold, palette,
                                    chr_colors = chr_colors,
                                    color_cycle = color_cycle,
                                    color_significant = color_significant,
                                    significant_color = significant_color,
                                    nonsignificant_color = nonsignificant_color,
                                    alpha = alpha)
  highlight_flag <- .gwas_plot_highlight(highlight, tab)
  out <- list(data = tab, chr = prep$chr, chr_labels = prep_chr_labels,
              chr_mid = as.numeric(prep$mids), threshold = threshold,
              cutline = threshold[1L], suggestive = suggestive,
              threshold_col = threshold_col, threshold_lty = threshold_lty,
              threshold_lwd = threshold_lwd, suggestive_col = suggestive_col,
              suggestive_lty = suggestive_lty, suggestive_lwd = suggestive_lwd,
              threshold_labels = threshold_labels,
              suggestive_labels = suggestive_labels,
              point_colors = point_colors, colors = point_colors,
              color_by = color_by, point_size = sizes$size,
              base_point_size = point_size, pve = sizes$pve, pve_max = sizes$max,
              point_size_by = point_size_by, highlight = highlight_flag,
              highlight_col = highlight_col, label_top = label_top,
              show_pve_legend = isTRUE(show_pve_legend), title = title,
              xlab = xlab, ylab = ylab, pch = pch, signal_cex = signal_cex,
              highlight_pch = highlight_pch, highlight_cex = highlight_cex,
              color_significant = isTRUE(color_significant),
              significant_color = significant_color,
              nonsignificant_color = nonsignificant_color, alpha = alpha,
              show_threshold_legend = isTRUE(show_threshold_legend),
              show_grid = isTRUE(show_grid), grid_color = grid_color,
              grid_lty = grid_lty, grid_lwd = grid_lwd, axis_cex = axis_cex,
              label_cex = label_cex, title_cex = title_cex,
              label_angle = as.numeric(label_angle[1L]),
              label_offset = label_offset, x_axis_las = as.numeric(x_axis_las[1L]),
              legend_cex = legend_cex, xlim = xlim, ylim = ylim,
              width_mode = prep$width_mode, chr_lengths = prep$chr_lengths,
              chr_spans = prep$spans, observed_chr_spans = prep$observed_spans,
              gap = prep$gap, gap_values = prep$gap_values,
              gap_mode = prep$gap_mode,
              call = match.call())
  class(out) <- c("manhattan_plot", "list")
  if (draw) print(out)
  invisible(out)
}

#' @export
print.manhattan_plot <- function(x, ...) {
  .draw_manhattan(x)
  invisible(x)
}

.prepare_qq <- function(tab, confidence) {
  if (!is.numeric(confidence) || length(confidence) != 1L ||
      !is.finite(confidence) || confidence <= 0 || confidence >= 1) {
    .stopf("confidence must be a single number between 0 and 1.")
  }
  keep <- is.finite(tab$p) & tab$p >= 0 & tab$p <= 1
  if (!any(keep)) .stopf("No P values in [0, 1] are available for a Q-Q plot.")
  tab <- tab[keep, , drop = FALSE]
  positive <- tab$p[tab$p > 0]
  floor_p <- if (length(positive)) max(.Machine$double.xmin, min(positive) / 10) else .Machine$double.xmin
  p <- pmax(pmin(tab$p, 1), floor_p)
  ord <- order(p, decreasing = TRUE)
  p_desc <- p[ord]
  n <- length(p_desc)
  i <- seq_len(n)
  expected_p <- (n - i + 0.5) / n
  expected <- -log10(expected_p)
  observed <- -log10(p_desc)
  alpha <- (1 - confidence) / 2
  order_i <- n - i + 1L
  lower_p <- stats::qbeta(alpha, order_i, i)
  upper_p <- stats::qbeta(1 - alpha, order_i, i)
  lower <- -log10(pmax(upper_p, .Machine$double.xmin))
  upper <- -log10(pmax(lower_p, .Machine$double.xmin))
  chi <- stats::qchisq(pmax(1 - p_desc, .Machine$double.xmin), df = 1)
  lambda <- stats::median(chi, na.rm = TRUE) / stats::qchisq(0.5, df = 1)
  list(data = tab, order = ord, expected = expected, observed = observed,
       lower = lower, upper = upper, lambda = lambda, n = n,
       floor_p = floor_p, confidence = confidence)
}

.draw_qq <- function(x) {
  xy <- c(x$expected, x$observed, x$upper)
  lim <- c(0, max(xy[is.finite(xy)], na.rm = TRUE) * 1.06)
  if (!is.finite(lim[2L]) || lim[2L] <= 0) lim[2L] <- 1
  graphics::plot(x$expected, x$observed, type = "n", xlim = lim, ylim = lim,
                 xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  graphics::polygon(c(x$expected, rev(x$expected)),
                    c(x$lower, rev(x$upper)),
                    col = grDevices::adjustcolor(x$band_col, alpha.f = 0.18),
                    border = NA)
  graphics::abline(0, 1, col = x$line_col, lwd = 1.1)
  graphics::points(x$expected, x$observed, pch = 16, cex = x$point_size,
                   col = x$point_col)
  if (any(x$highlight)) {
    jj <- which(x$highlight)
    graphics::points(x$expected[jj], x$observed[jj], pch = 21,
                     cex = x$point_size[jj] * 1.25, bg = x$highlight_col,
                     col = "white", lwd = 0.7)
  }
  graphics::axis(1, las = 1, cex.axis = 0.75)
  graphics::axis(2, las = 1, cex.axis = 0.75)
  graphics::mtext(x$xlab, side = 1, line = 1.55, cex = 0.82)
  graphics::mtext(x$ylab, side = 2, line = 2.45, cex = 0.84)
  if (isTRUE(x$show_lambda)) {
    graphics::legend("topleft", legend = sprintf("lambda = %.3f\nn = %d", x$lambda, x$n),
                     bty = "n", cex = 0.70, text.col = "#333333", inset = c(0.01, 0.01))
  }
  if (!is.null(x$title) && nzchar(x$title)) {
    graphics::mtext(x$title, side = 3, line = 0.65, cex = 1.04, font = 2)
  }
  graphics::box(col = "#666666")
  invisible(x)
}

#' Draw a publication-ready GWAS Q-Q plot
#'
#' @param gwas A GWAS table, a path accepted by [read_gwas()], or the result of
#'   [gwas_mlm()].
#' @param confidence Confidence level for the null order-statistic envelope.
#' @param point_col,band_col,line_col Point, confidence-band, and diagonal colors.
#' @param point_size Base point size.
#' @param highlight SNP IDs, `chr:position` strings, or an ID/position table.
#' @param highlight_col Color used for highlighted SNPs.
#' @param show_lambda Display genomic inflation lambda and the number of tests.
#' @param title,xlab,ylab Plot title and axis labels.
#' @param draw Draw immediately; otherwise return an object for later printing.
#' @return An object of class `qq_plot`.
#' @export
plot_qq <- function(gwas, confidence = 0.95, point_col = "#2C7FB8",
                    band_col = "#7B3294", line_col = "#333333", point_size = 0.62,
                    highlight = NULL, highlight_col = "#D73027",
                    show_lambda = TRUE, title = "GWAS Q-Q plot",
                    xlab = "Expected -log10(P)", ylab = "Observed -log10(P)",
                    draw = TRUE) {
  tab <- .gwas_plot_table(gwas)
  prep <- .prepare_qq(tab, confidence)
  highlight_flag <- .gwas_plot_highlight(highlight, prep$data)
  # `ord` indexes the retained table in decreasing P order.  Carry the flag
  # into the same order used by the plotted observed values.
  pve <- .gwas_point_sizes(prep$data, point_size, "PVE")
  out <- c(prep, list(point_col = point_col, band_col = band_col,
                      line_col = line_col, point_size = rep(point_size, prep$n),
                      highlight = highlight_flag[prep$order],
                      highlight_col = highlight_col, show_lambda = isTRUE(show_lambda),
                      title = title, xlab = xlab, ylab = ylab,
                      call = match.call()))
  class(out) <- c("qq_plot", "list")
  if (draw) print(out)
  invisible(out)
}

#' Save a Manhattan or Q-Q plot
#'
#' @param plot A `manhattan_plot` or `qq_plot` object.
#' @param file Output PDF, SVG, PNG, or TIFF file.
#' @param width,height Size in inches.
#' @param dpi Raster resolution.
#' @param bg Background color.
#' @return The normalized output path invisibly.
#' @export
save_gwas_plot <- function(plot, file, width = 7.2, height = 5.4,
                           dpi = 300, bg = "white") {
  if (!inherits(plot, c("manhattan_plot", "qq_plot"))) {
    .stopf("plot must be returned by plot_manhattan() or plot_qq().")
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(file))
  if (ext == "pdf") grDevices::pdf(file, width = width, height = height, onefile = TRUE, bg = bg)
  else if (ext == "svg") grDevices::svg(file, width = width, height = height, bg = bg)
  else if (ext == "png") grDevices::png(file, width = width, height = height,
                                         units = "in", res = dpi, bg = bg)
  else if (ext %in% c("tif", "tiff")) grDevices::tiff(file, width = width, height = height,
                                                        units = "in", res = dpi, bg = bg,
                                                        compression = "lzw")
  else .stopf("Unsupported plot extension '%s'. Use pdf, svg, png, tif, or tiff.", ext)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(normalizePath(file, mustWork = FALSE))
}

#' @export
print.qq_plot <- function(x, ...) {
  .draw_qq(x)
  invisible(x)
}
