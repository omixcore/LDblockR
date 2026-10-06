% Generated manually for LDblockR 0.0.1
\name{gwas_plots}
\alias{plot_manhattan}
\alias{plot_qq}
\alias{save_gwas_plot}
\alias{print.manhattan_plot}
\alias{print.qq_plot}
\title{Publication-ready Manhattan and Q-Q plots}
\description{
  Draws Manhattan and Q-Q diagnostics from a normalized GWAS table or the
  result of \code{gwas_mlm()}, and saves either plot in common publication
  formats without external plotting packages.
}
\usage{
plot_manhattan(gwas, cutline = NULL, suggestive = NULL,
  palette = c("#2C7FB8", "#7B3294"),
  color_by = c("chromosome", "significance", "pve"),
  point_size_by = c("PVE", "none"), point_size = 0.58,
  highlight = NULL, highlight_col = "#D73027", label_top = 0L,
  show_pve_legend = TRUE, title = "Manhattan plot", xlab = "Chromosome",
  ylab = expression(-log[10](P)), draw = TRUE, threshold = NULL,
  threshold_col = NULL, threshold_lty = 2, threshold_lwd = 1,
  threshold_color = NULL, suggestive_col = "#8C8C8C", suggestive_lty = 3,
  suggestive_lwd = 0.9, suggestive_color = NULL, chr_colors = NULL,
  color_cycle = TRUE, chr_labels = NULL, colors = NULL, gap_ratio = 0.05,
  gap_bp = NULL, gap_mode = c("fixed", "relative"), chr_lengths = NULL,
  width_mode = c("length", "observed", "equal"), ylim = NULL, xlim = NULL,
  pch = 16, point_cex = NULL,
  signal_cex = 1, color_significant = FALSE,
  significant_color = "#E64B35", nonsignificant_color = "#2C7FB8",
  alpha = 0.82, alpha_points = NULL, line_width = NULL,
  threshold_labels = NULL, suggestive_labels = NULL,
  show_threshold_legend = TRUE, show_grid = TRUE, grid_color = "#EEF2F5",
  grid_lty = 1, grid_lwd = 0.7, axis_cex = 0.72, label_cex = 0.58,
  title_cex = 1.04, label_angle = 35, label_offset = 0.025,
  x_axis_las = 1, highlight_pch = 21, highlight_cex = 1.25,
  legend_cex = 0.62)

plot_qq(gwas, confidence = 0.95, point_col = "#2C7FB8",
  band_col = "#7B3294", line_col = "#333333", point_size = 0.62,
  highlight = NULL, highlight_col = "#D73027", show_lambda = TRUE,
  title = "GWAS Q-Q plot", xlab = "Expected -log10(P)",
  ylab = "Observed -log10(P)", draw = TRUE)

save_gwas_plot(plot, file, width = 7.2, height = 5.4, dpi = 300,
  bg = "white")
}
\arguments{
  \item{gwas}{A GWAS table, path accepted by \code{read_gwas()}, or a
    \code{gwas_mlm_result}.}
  \item{cutline}{A \code{-log10(P)} threshold. \code{NULL} uses
    \code{-log10(0.05/n)}; \code{Inf} suppresses the line. Kept for
    backward compatibility; use \code{threshold} for multiple lines.}
  \item{suggestive}{Optional vector of suggestive \code{-log10(P)} thresholds.}
  \item{palette}{Fallback chromosome colors, recycled when
    \code{color_cycle = TRUE}.}
  \item{color_by}{Manhattan colors: \code{"chromosome"},
    \code{"significance"}, or \code{"pve"}.}
  \item{point_size_by}{Use marker-level PVE for Manhattan point sizes or use
    \code{"none"}.}
  \item{point_size}{Base point size.}
  \item{highlight}{SNP IDs, \code{chr:position} strings, or a table with ID
    and/or chromosome-position columns.}
  \item{highlight_col}{Highlight color.}
  \item{label_top}{Number of top Manhattan SNPs to label.}
  \item{show_pve_legend}{Whether to draw the PVE size legend.}
  \item{threshold}{One or more \code{-log10(P)} thresholds; overrides
    \code{cutline}.}
  \item{threshold_col,threshold_color}{Threshold-line colors; the latter is
    a compatibility alias.}
  \item{threshold_lty,threshold_lwd}{Threshold-line types and widths.}
  \item{suggestive_col,suggestive_color}{Suggestive-line colors.}
  \item{suggestive_lty,suggestive_lwd}{Suggestive-line types and widths.}
  \item{chr_colors}{Named or unnamed per-chromosome colors. Named values are
    matched to chromosome labels.}
  \item{color_cycle}{Recycle \code{chr_colors} or \code{palette} over
    chromosomes; \code{FALSE} requires enough colors.}
  \item{chr_labels}{Optional named or positional chromosome-axis labels.}
  \item{colors}{Alias for \code{chr_colors}.}
  \item{gap_ratio}{Scaling factor used to derive an automatic gap when
    \code{gap_bp} is not supplied. With the default \code{gap_mode = "fixed"},
    the resulting single gap is used at every chromosome boundary.}
  \item{gap_bp}{Explicit fixed gap on the genomic-position scale between
    chromosomes, in the same coordinate unit as \code{pos} (for example,
    \code{5e6} for a 5-Mb gap).}
  \item{gap_mode}{Gap rule: \code{"fixed"} uses one identical gap at every
    chromosome boundary (the publication default); \code{"relative"} scales
    each boundary gap by the width of the preceding chromosome. Supplying
    \code{gap_bp} always selects the fixed rule.}
  \item{chr_lengths}{Optional named or positional physical chromosome lengths
    in bp. They determine relative chromosome widths.}
  \item{width_mode}{Width rule when \code{chr_lengths} is absent:
    \code{"length"} uses the largest observed coordinate, \code{"observed"}
    uses the marker span, and \code{"equal"} gives every chromosome the same
    width.}
  \item{ylim,xlim}{Optional two-element plotting ranges.}
  \item{pch}{Point plotting symbol.}
  \item{point_cex}{Alias for \code{point_size}.}
  \item{signal_cex}{Size multiplier for points at or above the first threshold.}
  \item{color_significant}{Whether to recolor significant points in
    chromosome-color mode.}
  \item{significant_color,nonsignificant_color}{Colors for significance mode.}
  \item{alpha,alpha_points}{Point transparency in \code{[0, 1]}; the latter
    is an alias.}
  \item{line_width}{Set both threshold and suggestive line widths.}
  \item{threshold_labels,suggestive_labels}{Legend labels for threshold lines.}
  \item{show_threshold_legend}{Whether to draw a line legend.}
  \item{show_grid,grid_color,grid_lty,grid_lwd}{Horizontal grid-line controls.}
  \item{axis_cex,label_cex,title_cex}{Font scaling controls.}
  \item{label_angle,label_offset}{Top-SNP label rotation and vertical offset.}
  \item{x_axis_las}{Orientation of chromosome-axis labels.}
  \item{highlight_pch,highlight_cex}{Highlight symbol and size multiplier.}
  \item{legend_cex}{Legend font scaling.}
  \item{confidence}{Q-Q order-statistic confidence level between 0 and 1.}
  \item{point_col,band_col,line_col}{Q-Q point, confidence-band, and diagonal
    colors.}
  \item{show_lambda}{Whether to show genomic inflation lambda and test count.}
  \item{title,xlab,ylab}{Plot title and axis labels.}
  \item{draw}{Whether to draw immediately.}
  \item{plot}{A \code{manhattan_plot} or \code{qq_plot} object.}
  \item{file}{Output PDF, SVG, PNG, or TIFF path.}
  \item{width,height}{Output dimensions in inches.}
  \item{dpi}{Raster resolution for PNG/TIFF.}
  \item{bg}{Background color.}
}
\value{
  \code{plot_manhattan()} returns an object of class \code{manhattan_plot};
  \code{plot_qq()} returns an object of class \code{qq_plot};
  \code{save_gwas_plot()} returns the normalized output path invisibly.
}
\details{
  Manhattan coordinates are concatenated by chromosome with a small visual
  gap while retaining the original base-pair positions in the returned table.
  Q-Q confidence limits use exact beta order-statistic quantiles. P values of
  zero are floored only for plotting and do not alter the supplied GWAS table.
}
\examples{
\dontrun{
paths <- example_data("tassel")
pheno <- read_tassel_phenotype(paths[["mdp_phenotype"]])
g <- gwas_mlm(paths[["mdp_genotype"]], pheno, trait = "EarHT",
              covariates = c("location", "Q1", "Q2", "Q3"),
              replicate = "expand", verbose = FALSE)
m <- plot_manhattan(g, color_by = "significance", draw = FALSE)
q <- plot_qq(g, draw = FALSE)
save_gwas_plot(m, "manhattan.pdf")
save_gwas_plot(q, "qq.pdf")
}
}
