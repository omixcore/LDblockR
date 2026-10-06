#' Draw an integrated regional GWAS and LD view
#'
#' This convenience wrapper configures [plot_ld()] as the integrated regional
#' view used in the LDblockR workflow figure. It combines the regional GWAS
#' signal, optional gene models, MAF, key-SNP markers, block boundaries, and an
#' inverted-triangle LD heatmap on one shared SNP coordinate system. The
#' lower tracks remain aligned to the exact variants in `ld$data$variants`.
#'
#' @param ld An `ld_result` object.
#' @param metric One or both LD metrics.
#' @param gwas Optional regional GWAS file or table.
#' @param genes Optional GFF3/GTF file or normalized annotation table.
#' @param blocks Optional block table from [detect_ld_blocks()] or a
#'   chromosome/start/end table.
#' @param tags Optional tag-SNP table from [select_tag_snps()].
#' @param lead Lead variant ID, `chr:position`, or position.
#' @param special Optional highlighted variants or named labels.
#' @param cutline Regional GWAS `-log10(P)` threshold.
#' @param title Figure title.
#' @param palette LD heatmap palette.
#' @param heatmap_colors Optional independent LD heatmap palette.
#' @param position_scale Use equally spaced SNP indices or physical positions.
#' @param render Heatmap rendering mode.
#' @param show_maf Add the MAF track.
#' @param show_ld_colors Color association points by LD with the lead SNP.
#' @param gwas_color_by Association point color rule.
#' @param show_pve_legend Add the PVE size legend when available.
#' @param show_snp_labels Add the dense SNP label track.
#' @param show_snp_connectors Draw dashed SNP-to-heatmap guides.
#' @param show_key_snps Add the lead/tag/block-boundary key track.
#' @param show_connectors Add a separate vertical GWAS connector track.
#' @param show_gene_names,show_region_label Gene-track display switches.
#' @param draw Draw immediately.
#' @param ... Additional arguments passed to [plot_ld()].
#' @return An object of class `ld_plot`.
#' @export
plot_ld_region <- function(ld, metric = "r2", gwas = NULL, genes = NULL,
                           blocks = NULL, tags = NULL, lead = NULL,
                           special = NULL, cutline = Inf,
                           title = "Regional GWAS and LD view",
                           palette = "publication", heatmap_colors = NULL,
                           position_scale = "index", render = "auto",
                           show_maf = TRUE, show_ld_colors = TRUE,
                           gwas_color_by = "ld", show_pve_legend = FALSE,
                           show_snp_labels = FALSE,
                           show_snp_connectors = TRUE,
                           show_key_snps = TRUE, show_connectors = FALSE,
                           show_gene_names = TRUE, show_region_label = TRUE,
                           draw = TRUE, ...) {
  plot_ld(
    ld = ld, metric = metric, gwas = gwas, genes = genes, blocks = blocks,
    tags = tags, lead = lead, special = special, cutline = cutline,
    title = title, palette = palette, heatmap_colors = heatmap_colors,
    position_scale = position_scale, render = render, show_maf = show_maf,
    show_ld_colors = show_ld_colors, gwas_color_by = gwas_color_by,
    show_pve_legend = show_pve_legend, show_snp_labels = show_snp_labels,
    show_snp_connectors = show_snp_connectors,
    show_key_snps = show_key_snps, show_connectors = show_connectors,
    show_gene_names = show_gene_names, show_region_label = show_region_label,
    draw = draw, ...
  )
}
