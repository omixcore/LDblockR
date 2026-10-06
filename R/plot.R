.ld_palette <- function(palette, n = 101L) {
  if (length(palette) == 1L) {
    palette <- switch(tolower(palette),
      classic = c("#FFFDF2", "#FFF2A8", "#FDBB55", "#E84A33", "#B40426"),
      publication = c("#FFFDF2", "#FFF2A8", "#FDBB55", "#E84A33", "#B40426"),
      red = c("#FFFFFF", "#FDB863", "#B2182B"),
      blue = c("#F7FBFF", "#6BAED6", "#08306B"),
      purple = c("#FFFFFF", "#C2A5CF", "#762A83"),
      viridis = c("#440154", "#21908C", "#FDE725"),
      .stopf("Unknown palette '%s'.", palette)
    )
  }
  if (length(palette) < 2L || anyNA(palette)) .stopf("palette needs at least two valid colors.")
  grDevices::colorRampPalette(palette)(n)
}

.metric_label <- function(metric) if (metric == "r2") expression(r^2) else expression(D * "'")

.variant_boundaries <- function(x) {
  n <- length(x)
  if (n == 1L) return(c(x - 0.5, x + 0.5))
  middle <- (x[-n] + x[-1L]) / 2
  c(x[1L] - (middle[1L] - x[1L]), middle,
    x[n] + (x[n] - middle[n - 1L]))
}

.plot_coordinates <- function(variants, position_scale) {
  pos <- variants$pos
  if (position_scale == "physical" && all(diff(pos) > 0)) {
    x <- pos
    label_pos <- pos
  } else {
    if (position_scale == "physical" && any(diff(pos) <= 0)) {
      .warnf("Duplicate or unsorted positions detected; using equally spaced variant coordinates for plotting.")
    }
    x <- seq_along(pos)
    label_pos <- pos
    position_scale <- "index"
  }
  list(x = x, boundaries = .variant_boundaries(x), position_scale = position_scale,
       genomic_positions = label_pos)
}

.map_position_to_x <- function(pos, coord) {
  if (coord$position_scale == "physical") return(as.numeric(pos))
  stats::approx(coord$genomic_positions, coord$x, xout = as.numeric(pos),
                rule = 2, ties = "ordered")$y
}

# Match association rows to the exact variants used by the LD matrix.  The
# marker ID is preferred, with a chromosome/position key as a robust fallback
# for association files that renamed markers during export.  Returning the
# LD index (rather than re-interpolated coordinates) is important: every
# upper-track point and lower-track annotation then shares one SNP coordinate.
.match_gwas_to_ld <- function(gwas, variants) {
  if (is.null(gwas) || !nrow(gwas)) {
    return(list(gwas = gwas, match_idx = integer(), dropped = integer()))
  }
  gwas <- as.data.frame(gwas, stringsAsFactors = FALSE)
  if (!all(c("chr", "pos", "logp") %in% names(gwas))) {
    gwas <- read_gwas(gwas)
  }
  v_id <- as.character(variants$id)
  v_key <- paste(.chr_key(variants$chr), as.numeric(variants$pos), sep = "\r")
  g_key <- paste(.chr_key(gwas$chr), as.numeric(gwas$pos), sep = "\r")
  idx_id <- match(as.character(gwas$id), v_id)
  idx_key <- match(g_key, v_key)
  idx <- idx_id
  replace <- is.na(idx)
  # Do not allow an ID collision at a different coordinate to misalign a
  # point; this is common when tables contain duplicate marker names.
  valid_id <- !is.na(idx) & is.finite(gwas$pos) &
    .same_chr(variants$chr[idx], gwas$chr) &
    as.numeric(variants$pos[idx]) == as.numeric(gwas$pos)
  replace <- !valid_id
  idx[replace] <- idx_key[replace]
  keep <- !is.na(idx) & is.finite(gwas$pos) & is.finite(gwas$logp)
  dropped <- which(!keep)
  list(gwas = gwas[keep, , drop = FALSE], match_idx = as.integer(idx[keep]),
       dropped = dropped)
}

.format_position_axis <- function(at) {
  span <- diff(range(at, finite = TRUE))
  if (span >= 2e6) list(labels = format(at / 1e6, digits = 4, trim = TRUE), unit = "Position (Mb)")
  else if (span >= 2e3) list(labels = format(at / 1e3, digits = 4, trim = TRUE), unit = "Position (kb)")
  else list(labels = format(at, scientific = FALSE, trim = TRUE), unit = "Position (bp)")
}

.position_ticks <- function(coord, variants, n = 5L) {
  pos <- suppressWarnings(as.numeric(variants$pos))
  keep <- is.finite(pos)
  if (!any(keep)) return(list(at = numeric(), labels = character(), unit = "Position"))
  pos <- pos[keep]
  if (coord$position_scale == "physical") {
    at <- pretty(range(pos), n = n)
    at <- at[at >= min(pos) & at <= max(pos)]
    x <- .map_position_to_x(at, coord)
  } else {
    k <- unique(round(seq(1L, length(pos), length.out = min(6L, length(pos)))))
    at <- pos[k]
    x <- coord$x[k]
  }
  axis_info <- .format_position_axis(at)
  list(at = x, labels = axis_info$labels, unit = axis_info$unit)
}

.color_values <- function(value, palette) {
  idx <- floor(pmin(1, pmax(0, value)) * (length(palette) - 1L)) + 1L
  out <- palette[idx]
  out[!is.finite(value)] <- "#00000000"
  out
}

.draw_vector_triangle <- function(mat, coord, palette, border, border_width, show_values) {
  n <- nrow(mat)
  b <- coord$boundaries
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      value <- mat[i, j]
      if (!is.finite(value)) next
      u <- c(b[i], b[i + 1L], b[i + 1L], b[i])
      v <- c(b[j], b[j], b[j + 1L], b[j + 1L])
      # The LD matrix is anchored at the common SNP baseline and expands
      # downwards, matching the conventional inverted-triangle LD diagram.
      graphics::polygon((u + v) / 2, (u - v) / 2,
                        col = .color_values(value, palette),
                        border = border, lwd = border_width)
      if (show_values) {
        graphics::text((coord$x[i] + coord$x[j]) / 2,
                       (coord$x[i] - coord$x[j]) / 2,
                       labels = sprintf("%.2f", value), cex = 0.48)
      }
    }
  }
}

.draw_raster_triangle <- function(mat, coord, palette, raster_width = 1400L,
                                  raster_height = 700L) {
  b <- coord$boundaries
  xmin <- b[1L]
  xmax <- tail(b, 1L)
  ymax <- (xmax - xmin) / 2
  nx <- max(300L, min(as.integer(raster_width), max(300L, 2L * nrow(mat))))
  ny <- max(150L, min(as.integer(raster_height), max(150L, nrow(mat))))
  xs <- seq(xmin, xmax, length.out = nx)
  # Raster rows are filled from the SNP baseline downwards.  `rasterImage()`
  # expects the first row at the top of the image, hence no row reversal below.
  ys <- seq(0, -ymax, length.out = ny)
  cols <- matrix("#00000000", nrow = ny, ncol = nx)
  for (r in seq_len(ny)) {
    depth <- abs(ys[r])
    u <- xs - depth
    v <- xs + depth
    i <- findInterval(u, b, all.inside = TRUE)
    j <- findInterval(v, b, all.inside = TRUE)
    valid <- u >= xmin & v <= xmax & j > i & i >= 1L & j <= nrow(mat)
    value <- rep(NA_real_, nx)
    if (any(valid)) value[valid] <- mat[cbind(i[valid], j[valid])]
    cols[r, ] <- .color_values(value, palette)
  }
  img <- grDevices::as.raster(cols)
  graphics::rasterImage(img, xmin, -ymax, xmax, 0, interpolate = FALSE)
}

.draw_gwas_track <- function(ld, gwas, coord, xlim, cutline, lead, special,
                             point_size, show_ld_colors,
                             color_by = "significance", show_pve_legend = TRUE,
                             show_connectors = TRUE) {
  v <- ld$data$variants
  matched <- .match_gwas_to_ld(gwas, v)
  gwas <- matched$gwas
  match_idx <- matched$match_idx
  if (!nrow(gwas)) .stopf("No GWAS rows match the variants in the LD region.")
  # Use the exact LD coordinates, never an interpolated position.  This is
  # what keeps the association point, connector, and heatmap SNP aligned.
  gwas_x <- coord$x[match_idx]
  if (is.null(lead)) {
    lead_gwas <- which.max(gwas$logp)
    lead_idx <- match_idx[lead_gwas]
  } else if (is.character(lead)) {
    lead_idx <- match(lead[1L], v$id)
    if (is.na(lead_idx) && grepl(":", lead[1L], fixed = TRUE)) {
      lead_pos <- suppressWarnings(as.numeric(tail(strsplit(lead[1L], ":", fixed = TRUE)[[1L]], 1L)))
      if (is.finite(lead_pos)) lead_idx <- which.min(abs(v$pos - lead_pos))
    }
    if (is.na(lead_idx)) {
      lead_gwas <- match(lead[1L], gwas$id)
      if (!is.na(lead_gwas)) lead_idx <- match_idx[lead_gwas]
    }
    if (is.na(lead_idx)) .stopf("Lead variant '%s' was not found in LD data.", lead[1L])
  } else {
    lead_idx <- which.min(abs(v$pos - as.numeric(lead[1L])))
  }
  lead_row <- match(lead_idx, match_idx)
  if (is.na(lead_row)) {
    lead_row <- which.max(gwas$logp)
    lead_idx <- match_idx[lead_row]
  }
  color_by <- match.arg(color_by, c("significance", "pve", "ld"))
  pve <- rep(NA_real_, nrow(gwas))
  pve_col <- .match_column(gwas, c("PVE", "pve", "var_explained"), FALSE)
  if (!is.null(pve_col)) pve <- suppressWarnings(as.numeric(gwas[[pve_col]]))
  colors <- rep("#2C7FB8", nrow(gwas))
  if (color_by == "ld" && show_ld_colors && !is.null(ld$r2)) {
    rv <- rep(NA_real_, nrow(gwas))
    ok <- !is.na(match_idx)
    rv[ok] <- ld$r2[lead_idx, match_idx[ok]]
    colors <- .color_values(rv, .ld_palette(c("#4575B4", "#91CF60", "#D73027")))
    colors[!is.finite(rv)] <- "#BDBDBD"
  } else if (color_by == "pve" && any(is.finite(pve))) {
    pmaxv <- max(pve, na.rm = TRUE)
    if (!is.finite(pmaxv) || pmaxv <= 0) pmaxv <- 1
    colors <- .color_values(pmin(1, pmax(0, pve / pmaxv)),
                            .ld_palette(c("#2C7FB8", "#91CF60", "#E31A1C")))
    colors[!is.finite(pve)] <- "#BDBDBD"
  } else if (color_by == "significance") {
    sig <- is.finite(cutline) & gwas$logp >= cutline
    colors <- ifelse(sig, "#E31A1C", "#2C7FB8")
  }
  pve_scale <- if (any(is.finite(pve))) max(pve, na.rm = TRUE) else NA_real_
  cex <- rep(point_size, nrow(gwas))
  if (is.finite(pve_scale) && pve_scale > 0) {
    cex <- point_size * (0.78 + 1.22 * sqrt(pmin(1, pmax(0, pve / pve_scale))))
    cex[!is.finite(cex)] <- point_size
  }
  range_values <- c(gwas$logp, if (is.finite(cutline)) cutline else numeric())
  ylim <- c(0, max(range_values, na.rm = TRUE) * 1.12)
  if (!is.finite(ylim[2L]) || ylim[2L] <= 0) ylim[2L] <- 1
  graphics::plot(gwas_x, gwas$logp, type = "n", xlim = xlim, ylim = ylim,
                 xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  graphics::axis(2, las = 1, cex.axis = 0.75)
  graphics::mtext(expression(-log[10](P)), side = 2, line = 2.5, cex = 0.8)
  graphics::points(gwas_x, gwas$logp, pch = 16, cex = cex, col = colors)
  if (isTRUE(show_connectors)) {
    sig <- is.finite(cutline) & gwas$logp >= cutline
    stem_col <- ifelse(sig, "#D73027", "#333333")
    graphics::segments(gwas_x, gwas$logp, gwas_x, 0,
                       col = grDevices::adjustcolor(stem_col, alpha.f = 0.86),
                       lwd = ifelse(sig, 1.8, 1.35))
    graphics::points(gwas_x, gwas$logp, pch = 16, cex = cex, col = colors)
  }
  if (is.finite(cutline)) graphics::abline(h = cutline, lty = 2, col = "#D73027", lwd = 1)
  graphics::points(gwas_x[lead_row], gwas$logp[lead_row], pch = 23, bg = "#7B2CBF",
                   col = "white", cex = point_size * 1.8, lwd = 0.8)
  if (!is.null(special) && nrow(special)) {
    sp <- special[special$pos >= min(v$pos) & special$pos <= max(v$pos), , drop = FALSE]
    if (nrow(sp)) {
      for (i in seq_len(nrow(sp))) {
        row <- which.min(abs(gwas$pos - sp$pos[i]))
        graphics::segments(gwas_x[row], gwas$logp[row], gwas_x[row],
                           min(ylim[2], gwas$logp[row] + 0.07 * diff(ylim)), col = "#333333")
        graphics::text(gwas_x[row], min(ylim[2], gwas$logp[row] + 0.09 * diff(ylim)),
                       labels = sp$label[i], cex = 0.65, srt = 35, adj = c(0, 0.5))
      }
    }
  }
  if (show_pve_legend && is.finite(pve_scale) && pve_scale > 0) {
    lv <- c(0, pve_scale / 2, pve_scale)
    graphics::legend("topright", legend = formatC(lv, digits = 2, format = "f"),
                     title = "PVE", pch = 16,
                     pt.cex = point_size * (0.78 + 1.22 * sqrt(lv / pve_scale)),
                     col = "#4D4D4D", bty = "n", cex = 0.60,
                     inset = c(0.01, 0.01))
  }
  ticks <- .position_ticks(coord, v)
  if (length(ticks$at)) {
    graphics::axis(1, at = ticks$at, labels = ticks$labels,
                   cex.axis = 0.66, col.axis = "#333333",
                   col.ticks = "#555555", tck = -0.018)
    graphics::mtext(ticks$unit, side = 1, line = 1.35, cex = 0.72)
  }
  graphics::box(col = "#666666")
  list(lead_idx = lead_idx, match_idx = match_idx, gwas = gwas,
       x = gwas_x, ylim = ylim)
}

.gene_primary_rows <- function(genes) {
  primary <- genes[genes$type == "gene", , drop = FALSE]
  if (!nrow(primary)) primary <- genes[genes$type %in% c("mRNA", "transcript"), , drop = FALSE]
  if (!nrow(primary) && nrow(genes)) {
    key <- ifelse(!is.na(genes$parent), genes$parent, genes$id)
    groups <- split(seq_len(nrow(genes)), key)
    primary <- do.call(rbind, lapply(names(groups), function(nm) {
      z <- genes[groups[[nm]], , drop = FALSE]
      data.frame(chr = z$chr[1L], source = z$source[1L], type = "gene",
                 start = min(z$start), end = max(z$end), score = NA_real_,
                 strand = z$strand[1L], phase = ".", id = nm, parent = NA_character_,
                 name = nm, attributes = "", stringsAsFactors = FALSE)
    }))
  }
  primary
}

.assign_gene_lanes <- function(genes, max_lanes = 6L) {
  genes <- genes[order(genes$start, genes$end), , drop = FALSE]
  ends <- rep(-Inf, max_lanes)
  lane <- integer(nrow(genes))
  for (i in seq_len(nrow(genes))) {
    available <- which(ends < genes$start[i])
    z <- if (length(available)) available[1L] else which.min(ends)
    lane[i] <- z
    ends[z] <- max(ends[z], genes$end[i])
  }
  genes$lane <- lane
  genes
}

.draw_gene_track <- function(genes, coord, xlim, max_lanes, show_gene_names,
                             gene_colors, show_region_label = TRUE,
                             gwas_info = NULL, cutline = Inf) {
  primary <- .assign_gene_lanes(.gene_primary_rows(genes), max_lanes)
  nlanes <- max(1L, max(primary$lane, na.rm = TRUE))
  # Leave enough headroom for the Gene Name guide and for gene labels.  The
  # previous limit placed the guide on the top clipping boundary, which made
  # it appear truncated when the plot was saved to PDF/PNG.
  graphics::plot(NA, NA, xlim = xlim, ylim = c(0.35, nlanes + 1.10),
                 xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  # Carry the same SNP-wise connectors through the gene-model panel.  This
  # keeps the upper association points visibly linked to the lower heatmap
  # even when a gene track is inserted between the two panels.
  if (!is.null(gwas_info) && length(gwas_info$match_idx)) {
    idx <- gwas_info$match_idx
    keep <- is.finite(idx) & idx >= 1L & idx <= length(coord$x)
    idx <- idx[keep]
    if (length(idx)) {
      gw <- gwas_info$gwas[keep, , drop = FALSE]
      sig <- is.finite(cutline) & gw$logp >= cutline
      col <- ifelse(sig, "#D73027", "#98A2B3")
      graphics::segments(coord$x[idx], 0.42, coord$x[idx], nlanes + 0.68,
                         col = grDevices::adjustcolor(col, alpha.f = 0.38),
                         lwd = ifelse(sig, 1.0, 0.55), lty = 3)
    }
  }
  if (!nrow(primary)) {
    graphics::text(mean(xlim), 0.7, "No gene annotation in region", col = "#777777", cex = 0.75)
    return(invisible(NULL))
  }
  # A compact locus ribbon mirrors the reference layout: adjacent colored
  # intervals provide a visual gene/feature ruler, with connector lines to the
  # detailed models below.  The exact feature coordinates remain in the gene
  # model rows, so this ribbon is purely a readable overview.
  strip_y <- nlanes + 0.48
  strip_h <- 0.13
  graphics::rect(xlim[1L], strip_y - strip_h, xlim[2L], strip_y + strip_h,
                 col = "#A6CEE3", border = "#777777", lwd = 0.45)
  strip_cols <- c("#F39C12", "#F4B6C2", "#F1F509", "#A6CEE3", "#FFB000")
  for (i in seq_len(nrow(primary))) {
    xs0 <- .map_position_to_x(primary$start[i], coord)
    xe0 <- .map_position_to_x(primary$end[i], coord)
    graphics::rect(xs0, strip_y - strip_h, xe0, strip_y + strip_h,
                   col = strip_cols[(i - 1L) %% length(strip_cols) + 1L], border = NA)
    mid0 <- (xs0 + xe0) / 2
    graphics::segments(mid0, strip_y - strip_h, mid0, primary$lane[i] + 0.19,
                       col = gene_colors["arrow"], lwd = 0.65)
  }
  graphics::text(xlim[1L] + 0.012 * diff(xlim), strip_y + strip_h + 0.16, "Gene Name",
                 adj = c(0, 0), cex = 0.62)
  for (i in seq_len(nrow(primary))) {
    g <- primary[i, ]
    y <- g$lane
    xs <- .map_position_to_x(g$start, coord)
    xe <- .map_position_to_x(g$end, coord)
    graphics::segments(xs, y, xe, y, col = gene_colors["intron"], lwd = 1.4)
    transcript_ids <- genes$id[genes$type %in% c("mRNA", "transcript") &
                                 !is.na(genes$parent) & genes$parent == g$id]
    related <- (!is.na(genes$parent) & genes$parent %in% c(g$id, g$parent, transcript_ids)) |
      (!is.na(genes$id) & genes$id == g$id)
    child <- genes[genes$start <= g$end & genes$end >= g$start & related, , drop = FALSE]
    child <- child[child$type %in% c("CDS", "exon", "five_prime_UTR", "three_prime_UTR", "UTR"), , drop = FALSE]
    if (!nrow(child)) {
      graphics::rect(xs, y - 0.12, xe, y + 0.12, col = gene_colors["gene"], border = NA)
    } else {
      for (q in seq_len(nrow(child))) {
        cxs <- .map_position_to_x(max(child$start[q], g$start), coord)
        cxe <- .map_position_to_x(min(child$end[q], g$end), coord)
        col <- if (child$type[q] == "CDS") gene_colors["CDS"] else gene_colors["UTR"]
        height <- if (child$type[q] == "CDS") 0.18 else 0.11
        graphics::rect(cxs, y - height, cxe, y + height, col = col, border = NA)
      }
    }
    arrow_x <- if (g$strand == "-") xs else xe
    direction <- if (g$strand == "-") -1 else 1
    aw <- 0.012 * diff(xlim)
    graphics::polygon(c(arrow_x, arrow_x - direction * aw, arrow_x - direction * aw),
                      c(y, y + 0.16, y - 0.16), col = gene_colors["arrow"], border = NA)
    if (show_gene_names) {
      label <- if (!is.na(g$name) && nzchar(g$name)) g$name else g$id
      graphics::text((xs + xe) / 2, y + 0.30, label, cex = 0.62, adj = c(0.5, 0))
    }
  }
  graphics::mtext("Genes", side = 2, line = 1.2, cex = 0.75)
  if (show_region_label) {
    pos <- range(coord$genomic_positions, finite = TRUE)
    span <- diff(pos)
    if (!is.finite(span) || span <= 0) span <- 1
    unit <- if (span >= 1e6) "Mb" else if (span >= 1e3) "kb" else "bp"
    div <- if (unit == "Mb") 1e6 else if (unit == "kb") 1e3 else 1
    region_name <- if (nrow(genes) && "chr" %in% names(genes)) as.character(genes$chr[1L]) else "region"
    graphics::mtext(sprintf("Start:%s    %s Region:%.2f %s    End:%s",
                            format(pos[1L] / div, digits = 7, trim = TRUE),
                            region_name, span / div, unit,
                            format(pos[2L] / div, digits = 7, trim = TRUE)),
                    side = 1, line = 0.15, cex = 0.61)
  }
  invisible(primary)
}

.normalize_special <- function(special, variants) {
  if (is.null(special)) return(NULL)
  if (is.character(special) && is.null(names(special))) {
    idx <- match(special, variants$id)
    if (anyNA(idx)) .stopf("Special variants not found: %s", paste(special[is.na(idx)], collapse = ", "))
    return(data.frame(id = special, pos = variants$pos[idx], label = special, stringsAsFactors = FALSE))
  }
  if (is.character(special) && !is.null(names(special))) {
    idx <- match(names(special), variants$id)
    if (anyNA(idx)) .stopf("Special variants not found: %s", paste(names(special)[is.na(idx)], collapse = ", "))
    return(data.frame(id = names(special), pos = variants$pos[idx], label = unname(special), stringsAsFactors = FALSE))
  }
  sp <- as.data.frame(special, stringsAsFactors = FALSE)
  id_col <- .match_column(sp, c("id", "snp", "marker"), FALSE)
  pos_col <- .match_column(sp, c("pos", "position", "bp", "site"), FALSE)
  label_col <- .match_column(sp, c("label", "name", "annotation"), FALSE)
  if (is.null(pos_col) && is.null(id_col)) .stopf("special needs an ID or position column.")
  pos <- if (!is.null(pos_col)) as.numeric(sp[[pos_col]]) else variants$pos[match(sp[[id_col]], variants$id)]
  id <- if (!is.null(id_col)) as.character(sp[[id_col]]) else paste0(variants$chr[1L], ":", pos)
  label <- if (!is.null(label_col)) as.character(sp[[label_col]]) else id
  data.frame(id = id, pos = pos, label = label, stringsAsFactors = FALSE)
}

.draw_color_key <- function(x, y, width, height, palette, label) {
  n <- length(palette)
  edges <- seq(x, x + width, length.out = n + 1L)
  for (i in seq_len(n)) {
    graphics::rect(edges[i], y, edges[i + 1L], y + height,
                   col = palette[i], border = NA, xpd = NA)
  }
  graphics::rect(x, y, x + width, y + height, border = "#555555", lwd = 0.5, xpd = NA)
  graphics::text(x + width / 2, y + height * 1.85, label,
                 cex = 0.62, font = 2, xpd = NA)
  graphics::text(c(x, x + width / 2, x + width),
                 y - height * 0.65, c("0", "0.5", "1"), cex = 0.58, xpd = NA)
}

.plot_block_spans <- function(blocks, coord) {
  empty <- data.frame(block_index = integer(), label = character(),
                      left = numeric(), right = numeric(), col = character(),
                      stringsAsFactors = FALSE)
  if (is.null(blocks) || !nrow(blocks)) return(empty)
  has_indices <- all(c("start_index", "end_index") %in% names(blocks))
  has_coordinates <- all(c("start", "end") %in% names(blocks))
  if (!has_indices && !has_coordinates) {
    .warnf("Ignoring blocks: expected start_index/end_index or start/end columns.")
    return(empty)
  }
  labels <- if ("block" %in% names(blocks)) as.character(blocks$block) else
    paste0("Block", seq_len(nrow(blocks)))
  block_cols <- rep(c("#2166AC", "#B2182B", "#4D9221", "#762A83"),
                    length.out = nrow(blocks))
  rows <- list()
  k <- 0L
  for (i in seq_len(nrow(blocks))) {
    if (has_indices) {
      si <- suppressWarnings(as.integer(blocks$start_index[i]))
      ei <- suppressWarnings(as.integer(blocks$end_index[i]))
      if (any(!is.finite(c(si, ei))) || si < 1L || ei < si ||
          ei > length(coord$x)) next
    } else {
      start <- suppressWarnings(as.numeric(blocks$start[i]))
      end <- suppressWarnings(as.numeric(blocks$end[i]))
      if (any(!is.finite(c(start, end)))) next
      si <- which.min(abs(coord$genomic_positions - start))
      ei <- which.min(abs(coord$genomic_positions - end))
      if (si > ei) {
        z <- si
        si <- ei
        ei <- z
      }
    }
    left <- coord$boundaries[si]
    right <- coord$boundaries[ei + 1L]
    if (any(!is.finite(c(left, right))) || right <= left) next
    k <- k + 1L
    rows[[k]] <- data.frame(block_index = i, label = labels[i],
                            left = left, right = right, col = block_cols[i],
                            stringsAsFactors = FALSE)
  }
  if (!k) return(empty)
  do.call(rbind, rows)
}

.draw_heatmap_block_labels <- function(block_spans, xlim, ymax) {
  if (is.null(block_spans) || !nrow(block_spans)) return(invisible(NULL))
  for (i in seq_len(nrow(block_spans))) {
    left <- block_spans$left[i]
    right <- block_spans$right[i]
    mid <- (left + right) / 2
    depth <- min(ymax, (right - left) / 2)
    col <- block_spans$col[i]
    graphics::polygon(c(left, mid, right), c(0, -depth, 0),
                      border = col, lwd = 1.8, xpd = NA)
    label <- block_spans$label[i]
    if (!is.na(label) && nzchar(label)) {
      label_y <- -0.34 * depth
      label_cex <- 0.62
      label_w <- graphics::strwidth(label, cex = label_cex, font = 2)
      label_h <- graphics::strheight(label, cex = label_cex, font = 2)
      pad <- max(0.015 * diff(xlim), 0.35 * label_h)
      graphics::rect(mid - label_w / 2 - pad, label_y - label_h / 2 - pad / 2,
                     mid + label_w / 2 + pad, label_y + label_h / 2 + pad / 2,
                     col = grDevices::adjustcolor("white", alpha.f = 0.72),
                     border = NA, xpd = NA)
      graphics::text(mid, label_y, label, col = col, cex = label_cex,
                     font = 2, xpd = NA)
    }
  }
  invisible(NULL)
}

.snp_plot_labels <- function(variants, label_type = "id_position") {
  ids <- as.character(variants$id)
  positions <- format(variants$pos, scientific = FALSE, trim = TRUE, big.mark = ",")
  if (label_type == "id") return(ids)
  if (label_type == "position") return(positions)
  paste(ids, positions, sep = "\n")
}

.draw_snp_label_panel <- function(variants, coord, xlim,
                                  label_type = "id_position",
                                  max_labels = 120L,
                                  show_connectors = TRUE) {
  graphics::plot(NA, NA, xlim = xlim, ylim = c(0, 1),
                 xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  marker_anchor <- list(
    x_ndc = graphics::grconvertX(coord$x, from = "user", to = "ndc"),
    y_ndc = graphics::grconvertY(0.08, from = "user", to = "ndc")
  )
  n <- nrow(variants)
  keep <- seq_len(n)
  if (is.finite(max_labels)) max_labels <- max(1L, floor(max_labels))
  if (is.finite(max_labels) && n > max_labels) {
    keep <- unique(round(seq(1L, n, length.out = max_labels)))
    .warnf("Showing SNP labels for %d of %d markers; increase max_snp_labels to label every SNP.",
           length(keep), n)
  }
  labels <- .snp_plot_labels(variants, label_type)[keep]
  # Keep a faint vertical guide for every marker.  The SNP, Block, and
  # heatmap panels use the same x coordinates, so these guides join a marker
  # to the corresponding LD block even when the label track is thinned.
  if (isTRUE(show_connectors)) {
    graphics::segments(coord$x, 0, coord$x, 1,
                       col = grDevices::adjustcolor("#2CA25F", alpha.f = 0.52),
                       lwd = 0.55, lty = 2, xpd = NA)
  }
  graphics::segments(coord$x[1L], 0.08, coord$x[length(coord$x)], 0.08,
                     col = "#555555", lwd = 0.65)
  graphics::segments(coord$x[keep], 0.08, coord$x[keep], 0.18,
                     col = "#777777", lwd = 0.55)
  graphics::points(coord$x[keep], rep(0.08, length(keep)), pch = 25,
                   bg = "#2C7FB8", col = "white", cex = 0.52, xpd = NA)
  graphics::text(coord$x[keep], rep(0.93, length(keep)), labels = labels,
                 srt = 90, adj = c(1, 0.5), cex = 0.40,
                 col = "#333333", xpd = NA)
  graphics::mtext("SNPs", side = 2, line = 1.2, cex = 0.72)
  invisible(marker_anchor)
}

.plot_variant_index <- function(value, variants) {
  if (is.null(value) || !length(value)) return(NA_integer_)
  value <- value[1L]
  idx <- match(as.character(value), as.character(variants$id))
  if (!is.na(idx)) return(as.integer(idx))
  if (is.character(value) && grepl(":", value, fixed = TRUE)) {
    z <- strsplit(value, ":", fixed = TRUE)[[1L]]
    pos <- suppressWarnings(as.numeric(tail(z, 1L)))
    chr <- paste(z[-length(z)], collapse = ":")
    if (is.finite(pos)) {
      exact <- which(as.character(variants$chr) == chr & variants$pos == pos)
      if (length(exact)) return(as.integer(exact[1L]))
      return(as.integer(which.min(abs(variants$pos - pos))))
    }
  }
  pos <- suppressWarnings(as.numeric(value))
  if (is.finite(pos)) return(as.integer(which.min(abs(variants$pos - pos))))
  NA_integer_
}

.block_boundary_indices <- function(blocks, variants, coord) {
  if (is.null(blocks) || !nrow(blocks)) return(integer())
  out <- integer()
  if (all(c("start_index", "end_index") %in% names(blocks))) {
    out <- c(suppressWarnings(as.integer(blocks$start_index)),
             suppressWarnings(as.integer(blocks$end_index)))
  } else if (all(c("start", "end") %in% names(blocks))) {
    out <- c(vapply(blocks$start, function(z) {
      z <- suppressWarnings(as.numeric(z))
      if (!is.finite(z)) NA_integer_ else which.min(abs(variants$pos - z))
    }, integer(1L)),
    vapply(blocks$end, function(z) {
      z <- suppressWarnings(as.numeric(z))
      if (!is.finite(z)) NA_integer_ else which.min(abs(variants$pos - z))
    }, integer(1L)))
  }
  out <- unique(out[is.finite(out) & out >= 1L & out <= nrow(variants)])
  as.integer(out)
}

.draw_block_panel <- function(blocks, coord, xlim, gwas_info = NULL,
                              cutline = Inf, show_connectors = TRUE,
                              block_spans = NULL,
                              show_snp_connectors = TRUE,
                              tags = NULL, lead_idx = NULL, special = NULL,
                              variants = NULL, show_key_snps = TRUE,
                              key_snp_colors = c(lead = "#7B2CBF",
                                                 tag = "#F28E2B",
                                                 block = "#2166AC",
                                                 special = "#D73027")) {
  graphics::plot(NA, NA, xlim = xlim, ylim = c(0, 1),
                 xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  marker_anchor <- list(
    x_ndc = graphics::grconvertX(coord$x, from = "user", to = "ndc"),
    y_ndc = graphics::grconvertY(0, from = "user", to = "ndc")
  )
  if (isTRUE(show_snp_connectors)) {
    graphics::segments(coord$x, 0, coord$x, 1,
                       col = grDevices::adjustcolor("#2CA25F", alpha.f = 0.52),
                       lwd = 0.55, lty = 2, xpd = NA)
  }
  if (isTRUE(show_connectors) && !isTRUE(show_snp_connectors) &&
      !is.null(gwas_info) && length(gwas_info$match_idx)) {
    idx <- gwas_info$match_idx
    keep <- is.finite(idx) & idx >= 1L & idx <= length(coord$x)
    idx <- idx[keep]
    if (length(idx)) {
      gw <- gwas_info$gwas[keep, , drop = FALSE]
      sig <- is.finite(cutline) & gw$logp >= cutline
      graphics::segments(coord$x[idx], 0.03, coord$x[idx], 0.96,
                         col = grDevices::adjustcolor(ifelse(sig, "#D73027", "#98A2B3"),
                                                      alpha.f = 0.28),
                         lwd = ifelse(sig, 0.9, 0.55), lty = 3)
    }
  }
  if (is.null(block_spans)) block_spans <- .plot_block_spans(blocks, coord)
  baseline <- 0.78
  if (nrow(block_spans)) {
    for (i in seq_len(nrow(block_spans))) {
      xs <- block_spans$left[i]
      xe <- block_spans$right[i]
      mid <- (xs + xe) / 2
      depth <- 0.50 + 0.10 * ((i - 1L) %% 3L)
      graphics::segments(xs, baseline, mid, depth, col = block_spans$col[i], lwd = 1.55)
      graphics::segments(mid, depth, xe, baseline, col = block_spans$col[i], lwd = 1.55)
      graphics::text(mid, max(0.38, depth - 0.10), block_spans$label[i],
                     cex = 0.58, col = block_spans$col[i], font = 2)
    }
  }
  if (isTRUE(show_key_snps) && !is.null(variants)) {
    n <- nrow(variants)
    graphics::segments(coord$x[1L], 0.16, coord$x[n], 0.16,
                       col = "#555555", lwd = 0.65)
    boundary_idx <- .block_boundary_indices(blocks, variants, coord)
    if (length(boundary_idx)) {
      graphics::points(coord$x[boundary_idx], rep(0.16, length(boundary_idx)),
                       pch = 25, bg = key_snp_colors[["block"]], col = "white",
                       cex = 0.82, xpd = NA)
    }
    tag_ids <- character()
    if (is.character(tags)) tag_ids <- as.character(tags)
    else if (!is.null(tags) && nrow(tags)) {
      tags <- as.data.frame(tags, stringsAsFactors = FALSE)
      tag_col <- .match_column(tags, c("tag", "id", "snp", "marker"), FALSE)
      if (!is.null(tag_col)) tag_ids <- as.character(tags[[tag_col]])
    }
    tag_idx <- match(tag_ids, as.character(variants$id))
    tag_idx <- unique(tag_idx[is.finite(tag_idx)])
    if (length(tag_idx)) {
      graphics::points(coord$x[tag_idx], rep(0.16, length(tag_idx)), pch = 16,
                       col = key_snp_colors[["tag"]], cex = 0.95, xpd = NA)
    }
    if (length(lead_idx) && is.finite(lead_idx) && lead_idx >= 1L && lead_idx <= n) {
      graphics::points(coord$x[lead_idx], 0.16, pch = 18,
                       col = key_snp_colors[["lead"]], cex = 1.15, xpd = NA)
    }
    if (!is.null(special) && nrow(special)) {
      special_idx <- vapply(special$pos, function(z) {
        z <- suppressWarnings(as.numeric(z))
        if (!is.finite(z)) NA_integer_ else which.min(abs(variants$pos - z))
      }, integer(1L))
      special_idx <- unique(special_idx[is.finite(special_idx)])
      if (length(special_idx)) {
        graphics::points(coord$x[special_idx], rep(0.16, length(special_idx)), pch = 8,
                         col = key_snp_colors[["special"]], cex = 1.15, xpd = NA)
      }
    }
    legend_labels <- character()
    legend_pch <- numeric()
    legend_col <- character()
    legend_bg <- character()
    if (length(lead_idx) && is.finite(lead_idx)) {
      legend_labels <- c(legend_labels, "Lead SNP")
      legend_pch <- c(legend_pch, 18)
      legend_col <- c(legend_col, key_snp_colors[["lead"]])
      legend_bg <- c(legend_bg, key_snp_colors[["lead"]])
    }
    if (length(tag_idx)) {
      legend_labels <- c(legend_labels, "Tag SNP")
      legend_pch <- c(legend_pch, 16)
      legend_col <- c(legend_col, key_snp_colors[["tag"]])
      legend_bg <- c(legend_bg, key_snp_colors[["tag"]])
    }
    if (length(boundary_idx)) {
      legend_labels <- c(legend_labels, "Block boundary")
      legend_pch <- c(legend_pch, 25)
      legend_col <- c(legend_col, key_snp_colors[["block"]])
      legend_bg <- c(legend_bg, key_snp_colors[["block"]])
    }
    if (!is.null(special) && nrow(special)) {
      legend_labels <- c(legend_labels, "Highlighted")
      legend_pch <- c(legend_pch, 8)
      legend_col <- c(legend_col, key_snp_colors[["special"]])
      legend_bg <- c(legend_bg, key_snp_colors[["special"]])
    }
    if (length(legend_labels)) {
      graphics::legend("topleft", legend = legend_labels, pch = legend_pch,
                       col = legend_col, pt.bg = legend_bg, bty = "n",
                       ncol = min(2L, length(legend_labels)), cex = 0.55,
                       horiz = FALSE, inset = c(0.01, 0.01),
                       title = "Key SNPs", title.cex = 0.58)
    }
  }
  if (!nrow(block_spans) && !isTRUE(show_key_snps)) return(invisible(marker_anchor))
  graphics::mtext(if (isTRUE(show_key_snps)) "Key SNPs / blocks" else "Blocks",
                  side = 2, line = 1.2, cex = 0.72)
  invisible(marker_anchor)
}

.draw_direct_snp_connectors <- function(marker_anchor, coord) {
  if (is.null(marker_anchor) || !length(marker_anchor$x_ndc) ||
      !length(marker_anchor$y_ndc)) {
    return(invisible(NULL))
  }
  n <- length(coord$x)
  if (length(marker_anchor$x_ndc) != n) return(invisible(NULL))
  # Convert both ends through device coordinates.  This is important when
  # the heatmap panel has different margins or a different aspect ratio from
  # the upper SNP/block panel; using raw `coord$x` for only one end caused
  # connectors to fan into the wrong heatmap positions.
  x_end_ndc <- graphics::grconvertX(coord$x, from = "user", to = "ndc")
  y_end_ndc <- graphics::grconvertY(rep(0, n), from = "user", to = "ndc")
  x_start <- graphics::grconvertX(marker_anchor$x_ndc, from = "ndc", to = "user")
  y_start <- graphics::grconvertY(rep(marker_anchor$y_ndc, n),
                                  from = "ndc", to = "user")
  x_end <- graphics::grconvertX(x_end_ndc, from = "ndc", to = "user")
  y_end <- graphics::grconvertY(y_end_ndc, from = "ndc", to = "user")
  keep <- is.finite(x_start) & is.finite(y_start) &
    is.finite(x_end) & is.finite(y_end)
  if (any(keep)) {
    graphics::segments(x_start[keep], y_start[keep], x_end[keep], y_end[keep],
                       col = grDevices::adjustcolor("#2CA25F", alpha.f = 0.34),
                       lwd = 0.60, lty = 2, xpd = NA)
  }
  invisible(NULL)
}

.draw_heatmap_panel <- function(ld, metric, coord, tags, special,
                                gwas_info, cutline, palette, render,
                                show_values, grid_color, grid_width,
                                raster_width, raster_height, block_spans = NULL,
                                marker_anchor = NULL,
                                show_snp_connectors = TRUE) {
  mat <- ld[[metric]]
  b <- coord$boundaries
  xlim <- range(b)
  ymax <- diff(xlim) / 2
  below <- 0.18 * diff(xlim)
  top <- 0.10 * diff(xlim)
  # Keep the SNP baseline at y = 0 and place the triangular LD matrix below
  # it, as in the reference regional-association diagrams.
  # Keep the heatmap on the same full-width x scale as every upper track.
  # `asp = 1` made the plotting region shrink to a centered narrow column
  # whenever the triangular matrix was taller than the panel; consequently
  # the heatmap was clipped visually and SNP connectors no longer aligned.
  # The triangle is drawn in the shared coordinate system below, so its
  # horizontal extent and all connector endpoints remain fully visible.
  graphics::plot(NA, NA, xlim = xlim, ylim = c(-ymax - below, top),
                 xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "",
                 asp = NA)
  # Connect the real SNP marker positions to the heatmap's actual top edge.
  # With square cells the heatmap can be narrower than the upper tracks, so
  # the endpoints are converted through device coordinates before drawing.
  if (isTRUE(show_snp_connectors)) {
    .draw_direct_snp_connectors(marker_anchor, coord)
  }
  if (render == "auto") render <- if (nrow(mat) <= 250L) "vector" else "raster"
  if (render == "vector") {
    .draw_vector_triangle(mat, coord, palette, grid_color, grid_width,
                          show_values && nrow(mat) <= 40L)
  } else {
    .draw_raster_triangle(mat, coord, palette, raster_width, raster_height)
  }
  # The baseline follows the first and last SNP coordinates rather than the
  # plotting margins, so the upper edge of the inverted triangle and every
  # connector terminate on the same genomic ruler.
  graphics::segments(coord$x[1L], 0, coord$x[length(coord$x)], 0,
                     col = "#333333", lwd = 0.9)
  # Draw the LDBlockShow-style site guides before the GWAS markers so that
  # significant GWAS markers remain visible above the green guides.
  if (isTRUE(show_snp_connectors)) {
    graphics::segments(coord$x, 0, coord$x, -0.055 * ymax,
                       col = grDevices::adjustcolor("#2CA25F", alpha.f = 0.52),
                       lwd = 0.55, lty = 2, xpd = NA)
  }
  # Repeat the GWAS marker coordinates at the heatmap baseline.  The short
  # stems continue the connector panel visually and make it obvious which
  # SNP in the triangular matrix corresponds to each association point.
  if (!is.null(gwas_info) && length(gwas_info$match_idx)) {
    idx <- gwas_info$match_idx
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= length(coord$x)]
    if (length(idx)) {
      gw <- gwas_info$gwas
      sig <- is.finite(cutline) & gw$logp >= cutline
      sig <- sig[seq_len(min(length(sig), length(gwas_info$match_idx)))]
      sig <- sig[match(idx, gwas_info$match_idx)]
      col <- ifelse(isTRUE(cutline == Inf), "#7B8794",
                    ifelse(sig, "#D73027", "#7B8794"))
      graphics::segments(coord$x[idx], 0, coord$x[idx], -0.065 * ymax,
                         col = col, lwd = 1.1, lty = 3, xpd = NA)
      graphics::points(coord$x[idx], rep(0, length(idx)), pch = 25,
                       bg = col, col = "white", cex = 0.75, xpd = NA)
    }
  }
  # Outline each within-block LD sub-triangle and label it in the heatmap.
  # These are block-specific guides, not a global triangle hull.
  .draw_heatmap_block_labels(block_spans, xlim, ymax)
  if (!is.null(tags) && nrow(tags)) {
    idx <- match(tags$tag, ld$data$variants$id)
    idx <- idx[!is.na(idx)]
    graphics::points(coord$x[idx], rep(0, length(idx)), pch = 25, bg = "#2B8CBE",
                     col = "white", cex = 0.8, xpd = NA)
  }
  if (!is.null(special) && nrow(special)) {
    sx <- .map_position_to_x(special$pos, coord)
    graphics::segments(sx, 0, sx, -0.045 * ymax, col = "#111111", lwd = 1)
    graphics::text(sx, -0.055 * ymax, special$label, srt = 90, adj = c(0, 0.5), cex = 0.55)
  }
  key_label <- if (metric == "r2") expression(R^2 ~ "color key") else
    expression(D * "'" ~ "color key")
  .draw_color_key(xlim[2L] - 0.28 * diff(xlim), -ymax - 0.13 * ymax,
                  0.24 * diff(xlim), 0.028 * ymax, palette,
                  key_label)
  graphics::mtext(if (metric == "r2") expression("LD heatmap (" * r^2 * ")") else expression("LD heatmap (" * D * "')"),
                  side = 2, line = 0.5, cex = 0.72)
}

.draw_connector_panel <- function(gwas_info, coord, xlim, cutline) {
  graphics::plot(NA, NA, xlim = xlim, ylim = c(0, 1), xaxs = "i", yaxs = "i",
                 axes = FALSE, xlab = "", ylab = "")
  if (!is.null(gwas_info) && length(gwas_info$match_idx)) {
    idx <- gwas_info$match_idx
    keep <- is.finite(idx) & idx >= 1L & idx <= length(coord$x)
    idx <- idx[keep]
    if (length(idx)) {
      gw <- gwas_info$gwas[keep, , drop = FALSE]
      sig <- is.finite(cutline) & gw$logp >= cutline
      col <- ifelse(sig, "#D73027", "#98A2B3")
      graphics::segments(coord$x[idx], 0.98, coord$x[idx], 0.02,
                         col = grDevices::adjustcolor(ifelse(sig, "#D73027", "#333333"),
                                                      alpha.f = 0.86),
                         lwd = ifelse(sig, 1.8, 1.35))
      graphics::points(coord$x[idx], rep(0.94, length(idx)), pch = 16,
                       col = ifelse(sig, "#D73027", "#333333"), cex = 0.45)
      graphics::points(coord$x[idx], rep(0.02, length(idx)), pch = 25,
                       bg = ifelse(sig, "#D73027", "#333333"),
                       col = "white", cex = 0.62)
    }
  }
  graphics::box(col = "#D0D5DD", lwd = 0.45)
}

.draw_ld_plot <- function(x) {
  ld <- x$ld
  v <- ld$data$variants
  coord <- .plot_coordinates(v, x$position_scale)
  xlim <- range(coord$boundaries)
  block_spans <- .plot_block_spans(x$blocks, coord)
  gwas_info <- x$gwas_match
  if (is.null(gwas_info) && !is.null(x$gwas) && nrow(x$gwas)) {
    gwas_info <- .match_gwas_to_ld(x$gwas, v)
  }
  lead_idx <- .plot_variant_index(x$lead, v)
  key_snp_colors <- x$key_snp_colors
  if (is.null(key_snp_colors)) {
    key_snp_colors <- c(lead = "#7B2CBF", tag = "#F28E2B",
                        block = "#2166AC", special = "#D73027")
  }
  panels <- character()
  heights <- numeric()
  if (!is.null(x$gwas) && nrow(x$gwas)) {
    panels <- c(panels, "gwas")
    heights <- c(heights, 2.25)
    if (isTRUE(x$show_connectors)) {
      panels <- c(panels, "connectors")
      heights <- c(heights, 0.38)
    }
  }
  if (!is.null(x$genes)) {
    panels <- c(panels, "genes")
      heights <- c(heights, 1.55)
  }
  if (isTRUE(x$show_maf)) {
    panels <- c(panels, "maf")
    heights <- c(heights, 0.9)
  }
  if (isTRUE(x$show_snp_labels)) {
    panels <- c(panels, "snp_labels")
    heights <- c(heights, 1.05)
  }
  tag_track <- !is.null(x$tags) &&
    (is.character(x$tags) && length(x$tags) > 0L ||
     !is.character(x$tags) && nrow(as.data.frame(x$tags)) > 0L)
  has_key_track <- isTRUE(x$show_key_snps) &&
    ((!is.null(x$gwas) && nrow(x$gwas)) || tag_track ||
     (!is.null(x$special) && nrow(x$special)) ||
     !is.na(lead_idx))
  if ((!is.null(x$blocks) && nrow(x$blocks)) || has_key_track) {
    panels <- c(panels, "blocks")
    heights <- c(heights, 0.85)
  }
  panels <- c(panels, paste0("heatmap:", x$metrics))
  heights <- c(heights, rep(3.55, length(x$metrics)))
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::layout(matrix(seq_along(panels), ncol = 1L), heights = heights)
  lead_idx <- .plot_variant_index(x$lead, v)
  marker_anchor <- NULL
  for (panel in panels) {
    if (panel == "gwas") {
      graphics::par(mar = c(2.45, 4.1, if (is.null(x$title)) 0.8 else 2.0, 1.0))
      gwas_draw <- .draw_gwas_track(ld, x$gwas, coord, xlim, x$cutline, x$lead,
                                    x$special, x$point_size, x$show_ld_colors,
                                    x$gwas_color_by, x$show_pve_legend,
                                    x$show_connectors)
      lead_idx <- gwas_draw$lead_idx
      if (!is.null(x$title)) graphics::mtext(x$title, side = 3, line = 0.7, cex = 1.05, font = 2)
    } else if (panel == "connectors") {
      graphics::par(mar = c(0.02, 4.1, 0.02, 1.0))
      .draw_connector_panel(gwas_info, coord, xlim, x$cutline)
    } else if (panel == "genes") {
      graphics::par(mar = c(1.10, 4.1, 0.30, 1.0))
      .draw_gene_track(x$genes, coord, xlim, x$max_gene_lanes,
                       x$show_gene_names, x$gene_colors, x$show_region_label,
                       gwas_info, x$cutline)
    } else if (panel == "maf") {
      graphics::par(mar = c(0.2, 4.1, 0.2, 1.0))
      graphics::plot(coord$x, v$maf, type = "h", lwd = 1.1, col = "#2C7FB8",
                     xlim = xlim, ylim = c(0, 0.5), xaxs = "i", yaxs = "i",
                     axes = FALSE, xlab = "", ylab = "")
      if (!is.null(gwas_info) && length(gwas_info$match_idx)) {
        idx <- gwas_info$match_idx
        keep <- is.finite(idx) & idx >= 1L & idx <= length(coord$x)
        if (any(keep)) graphics::segments(coord$x[idx[keep]], 0,
                                          coord$x[idx[keep]], 0.5,
                                          col = grDevices::adjustcolor("#98A2B3", alpha.f = 0.35),
                                          lty = 3, lwd = 0.7)
      }
      graphics::points(coord$x, v$maf, pch = 16, cex = 0.35, col = "#2C7FB8")
      graphics::axis(2, at = c(0, 0.25, 0.5), las = 1, cex.axis = 0.65)
      graphics::mtext("MAF", side = 2, line = 2.2, cex = 0.72)
      graphics::box(col = "#666666")
    } else if (panel == "snp_labels") {
      graphics::par(mar = c(0.02, 4.1, 0.02, 1.0))
      marker_anchor <- .draw_snp_label_panel(v, coord, xlim, x$snp_label_type,
                                              x$max_snp_labels, x$show_snp_connectors)
    } else if (panel == "blocks") {
      graphics::par(mar = c(0.05, 4.1, 0.08, 1.0))
      block_anchor <- .draw_block_panel(x$blocks, coord, xlim, gwas_info, x$cutline,
                                        x$show_connectors, block_spans,
                                        x$show_snp_connectors,
                                        tags = x$tags, lead_idx = lead_idx,
                                        special = x$special, variants = v,
                                        show_key_snps = isTRUE(x$show_key_snps),
                                        key_snp_colors = key_snp_colors)
      if (is.null(marker_anchor)) marker_anchor <- block_anchor
    } else {
      metric <- sub("^heatmap:", "", panel)
      graphics::par(mar = c(0.75, 4.1, 0.2, 1.0))
      .draw_heatmap_panel(ld, metric, coord, x$tags, x$special,
                          gwas_info, x$cutline, x$palette, x$render,
                          x$show_values, x$grid_color, x$grid_width,
                          x$raster_width, x$raster_height, block_spans,
                          marker_anchor,
                          x$show_snp_connectors)
      if (is.null(x$gwas) && is.null(x$title) == FALSE && metric == x$metrics[1L]) {
        graphics::mtext(x$title, side = 3, line = 0.4, cex = 1.05, font = 2)
      }
    }
  }
  invisible(lead_idx)
}

#' Draw an LD heatmap with aligned genomic tracks
#'
#' @param ld An `ld_result` object.
#' @param metric `"r2"`, `"dprime"`, `"both"`, or a vector of metrics.
#' @param gwas Optional GWAS file or data frame.
#' @param genes Optional GFF3 file or normalized annotation data frame.
#' @param blocks Optional block table. Blocks are drawn in a dedicated track
#' immediately above the LD heatmap; both `ld_blocks` output and `chr/start/end`
#' tables are accepted.
#' @param tags Optional tag-SNP table.
#' @param lead Lead variant ID or position; the top GWAS point is used by default.
#' @param special Highlighted variants, named labels, or a table.
#' @param cutline GWAS -log10(P) threshold.
#' @param title Plot title.
#' @param palette LD palette name or color vector.
#' @param heatmap_colors Optional LD heatmap palette name or color vector;
#'   when supplied, it overrides `palette` for the heatmap and color key.
#' @param position_scale Physical or equally spaced SNP coordinates. The
#'   default `"index"` uses equally spaced marker coordinates on the same
#'   full-width x scale as the upper tracks; use `"physical"` to scale SNPs
#'   by base-pair distance.
#' @param render Automatic, vector, or raster heatmap rendering.
#' @param show_values Print cell values for small vector heatmaps.
#' @param show_maf Include an allele-frequency track (off by default for the
#'   publication layout).
#' @param show_ld_colors Allow LD-based GWAS point colors when
#'   `gwas_color_by = "ld"`.
#' @param show_snp_labels Add an aligned SNP marker/label track above the
#'   Block track. The heatmap baseline also receives a small tick for every
#'   LD marker.
#' @param snp_label_type Label each SNP by its ID, physical position, or both.
#' @param max_snp_labels Maximum number of SNP labels. When fewer markers are
#'   present, every SNP is labeled; set to `Inf` to disable thinning.
#' @param show_snp_connectors Draw faint dashed guides from every SNP marker
#'   directly to the actual top edge of the LD heatmap.
#' @param show_key_snps Add a compact key-SNP track with lead-SNP, tag-SNP,
#'   block-boundary, and highlighted-variant markers.
#' @param key_snp_colors Named colors for `lead`, `tag`, `block`, and `special`
#'   markers in the key-SNP track.
#' @param gwas_color_by Color association points by significance, PVE, or LD.
#' @param show_pve_legend Add a point-size legend when the GWAS table has PVE.
#' @param show_connectors Draw aligned lines from each GWAS point to its LD
#'   marker at the top of the heatmap.
#' @param show_region_label Add Start/Region/End text to the gene track.
#' @param show_gene_names Label genes.
#' @param max_gene_lanes Maximum gene lanes.
#' @param point_size GWAS point size.
#' @param grid_color,grid_width Vector-cell border appearance.
#' @param gene_colors Named colors for gene components.
#' @param raster_width,raster_height Raster resolution for large heatmaps.
#' @param draw Draw immediately.
#' @return An object of class `ld_plot`.
#' @export
plot_ld <- function(ld, metric = "r2",
                    gwas = NULL, genes = NULL, blocks = NULL, tags = NULL,
                    lead = NULL, special = NULL, cutline = Inf, title = NULL,
                    palette = "classic", position_scale = c("index", "physical"),
                    render = c("auto", "vector", "raster"), show_values = FALSE,
                    show_maf = FALSE, show_ld_colors = FALSE,
  show_snp_labels = TRUE,
  snp_label_type = c("id_position", "id", "position"),
  max_snp_labels = 120L,
  show_snp_connectors = TRUE,
  show_key_snps = TRUE,
  key_snp_colors = c(lead = "#7B2CBF", tag = "#F28E2B",
                    block = "#2166AC", special = "#D73027"),
  gwas_color_by = c("significance", "pve", "ld"),
                    show_pve_legend = TRUE, show_connectors = TRUE,
                    show_region_label = TRUE,
                    show_gene_names = TRUE, max_gene_lanes = 6L,
                    point_size = 0.65, grid_color = "white", grid_width = 0.40,
                    gene_colors = c(CDS = "#D95F02", UTR = "#7570B3",
                                    intron = "#80B1D3", gene = "#1B9E77", arrow = "#333333"),
                    raster_width = 1400L, raster_height = 700L,
                    draw = TRUE, heatmap_colors = NULL) {
  if (!inherits(ld, "ld_result")) .stopf("ld must be an ld_result object.")
  position_scale <- match.arg(position_scale)
  render <- match.arg(render)
  snp_label_type <- match.arg(snp_label_type)
  gwas_color_by <- match.arg(gwas_color_by)
  max_snp_labels <- as.numeric(max_snp_labels)[1L]
  if (is.na(max_snp_labels) || max_snp_labels < 1) {
    .stopf("max_snp_labels must be a positive number or Inf.")
  }
  if (is.finite(max_snp_labels)) max_snp_labels <- max(1L, floor(max_snp_labels))
  key_snp_colors <- as.character(key_snp_colors)
  if (is.null(names(key_snp_colors)) && length(key_snp_colors) == 4L) {
    names(key_snp_colors) <- c("lead", "tag", "block", "special")
  }
  required_key_colors <- c("lead", "tag", "block", "special")
  if (is.null(names(key_snp_colors)) ||
      any(!required_key_colors %in% names(key_snp_colors)) ||
      anyNA(key_snp_colors[required_key_colors])) {
    .stopf("key_snp_colors must contain named lead, tag, block, and special colors.")
  }
  key_snp_colors <- key_snp_colors[required_key_colors]
  metric <- as.character(metric)
  if (!length(metric) || any(!metric %in% c("r2", "dprime", "both"))) {
    .stopf("metric must be r2, dprime, or both.")
  }
  if ("both" %in% metric) metric <- c("r2", "dprime") else metric <- unique(metric)
  for (z in metric) if (is.null(ld[[z]])) .stopf("%s was not calculated.", z)
  heatmap_palette_input <- if (is.null(heatmap_colors)) palette else heatmap_colors
  heatmap_palette <- .ld_palette(heatmap_palette_input)
  region <- paste0(ld$data$variants$chr[1L], ":", min(ld$data$variants$pos), "-", max(ld$data$variants$pos))
  if (is.character(gwas) && length(gwas) == 1L) gwas <- read_gwas(gwas, region = region)
  else if (!is.null(gwas)) gwas <- read_gwas(gwas, region = region)
  gwas_match <- NULL
  if (!is.null(gwas)) {
    gwas_match <- .match_gwas_to_ld(gwas, ld$data$variants)
    if (!nrow(gwas_match$gwas)) {
      .stopf("No GWAS rows match the variants in the LD region. Supply the regional GWAS table produced by gwas_ld_region().")
    }
    if (length(gwas_match$dropped)) {
      .warnf("Dropped %d GWAS rows that are not present in the LD matrix; plotting only exact SNP matches.",
             length(gwas_match$dropped))
    }
    gwas <- gwas_match$gwas
  }
  if (is.character(genes) && length(genes) == 1L) genes <- read_gff3(genes, region = region)
  else if (!is.null(genes)) genes <- as.data.frame(genes, stringsAsFactors = FALSE)
  if (is.character(tags)) tags <- data.frame(tag = as.character(tags), stringsAsFactors = FALSE)
  else if (!is.null(tags)) tags <- as.data.frame(tags, stringsAsFactors = FALSE)
  special <- .normalize_special(special, ld$data$variants)
  p <- list(
    ld = ld, metrics = metric, gwas = gwas, genes = genes, blocks = blocks,
    tags = tags, lead = lead, special = special, cutline = cutline, title = title,
    gwas_match = gwas_match, show_connectors = isTRUE(show_connectors),
    palette = heatmap_palette, heatmap_colors = heatmap_palette,
    position_scale = position_scale, render = render,
    show_values = show_values, show_maf = show_maf, show_ld_colors = show_ld_colors,
    show_snp_labels = isTRUE(show_snp_labels), snp_label_type = snp_label_type,
    max_snp_labels = max_snp_labels,
    show_snp_connectors = isTRUE(show_snp_connectors),
    show_key_snps = isTRUE(show_key_snps), key_snp_colors = key_snp_colors,
    gwas_color_by = gwas_color_by, show_pve_legend = show_pve_legend,
    show_region_label = show_region_label,
    show_gene_names = show_gene_names, max_gene_lanes = as.integer(max_gene_lanes),
    point_size = point_size, grid_color = grid_color, grid_width = grid_width,
    gene_colors = gene_colors, raster_width = raster_width, raster_height = raster_height
  )
  class(p) <- "ld_plot"
  if (draw) print(p)
  invisible(p)
}

#' @export
print.ld_plot <- function(x, ...) {
  .draw_ld_plot(x)
  invisible(x)
}

#' Save an LD plot
#'
#' @param plot An `ld_plot` object.
#' @param file Output PDF, SVG, PNG, or TIFF file.
#' @param width,height Size in inches.
#' @param dpi Raster resolution.
#' @param bg Background color.
#' @return The normalized output path invisibly.
#' @export
save_ld_plot <- function(plot, file, width = 11, height = 8.5, dpi = 300, bg = "white") {
  if (!inherits(plot, "ld_plot")) .stopf("plot must be returned by plot_ld().")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(file))
  if (ext == "pdf") grDevices::pdf(file, width = width, height = height, onefile = TRUE, bg = bg)
  else if (ext == "svg") grDevices::svg(file, width = width, height = height, bg = bg)
  else if (ext == "png") grDevices::png(file, width = width, height = height, units = "in", res = dpi, bg = bg)
  else if (ext %in% c("tif", "tiff")) grDevices::tiff(file, width = width, height = height, units = "in", res = dpi, bg = bg, compression = "lzw")
  else .stopf("Unsupported plot extension '%s'. Use pdf, svg, png, tif, or tiff.", ext)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(normalizePath(file, mustWork = FALSE))
}
